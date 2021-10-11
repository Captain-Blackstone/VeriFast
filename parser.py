import xml.etree.ElementTree as ET
from collections import defaultdict
import collections
from itertools import chain

import json
import unicodedata
import nltk

# nltk.download('punkt')


class OrderedDefaultdict(collections.OrderedDict):
    """ A defaultdict with OrderedDict as its base class. """

    def __init__(self, default_factory=None, *args, **kwargs):
        if not (default_factory is None or callable(default_factory)):
            raise TypeError('first argument must be callable or None')
        super(OrderedDefaultdict, self).__init__(*args, **kwargs)
        self.default_factory = default_factory  # called by __missing__()

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key, )
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):  # Optional, for pickle support.
        args = (self.default_factory,) if self.default_factory else tuple()
        return self.__class__, args, None, None, iter(self.items())

    def __repr__(self):  # Optional.
        return '%s(%r, %r)' % (self.__class__.__name__, self.default_factory, self.items())


class Parser:
    def __init__(self):
        """
        This is a Parser module. It is used to parse papers downloaded
        from different sources into the unified format.
        """
        pass

    def _work_with_section(self, section, paragraphs, tag=""):
        children = [section[x] for x in range(len(section))]
        for child in children:
            if child.tag in ("title", "label"):

                text = child.text
                if text is None:
                    if len(list(child.getchildren())) > 0:
                        text = child[0].text
                if text is not None:
                    if not tag:
                        tag = text
                    else:
                        tag += "|" + text

            if child.tag == "p":
                specital_paragraph = False
                paragraph_children = child.getchildren()
                for p_child in paragraph_children:
                    if p_child.tag in ("fig", "table-wrap"):
                        paragraphs = self._work_with_section(p_child, paragraphs, tag)
                        specital_paragraph = True
                if not specital_paragraph:
                    paragraph = "".join(child.itertext()).strip()
                else:
                    parts = ([child.text] +
                             list(chain(*([c.text, c.tail] for c in child.getchildren() if
                                          c.tag not in ("fig", "table-wrap")))) +
                             [child.tail])
                    paragraph = ''.join(filter(None, parts)).strip()
                if tag:
                    paragraphs[tag].append(paragraph)
                else:
                    paragraphs["Introduction"].append(paragraph)
            else:
                paragraphs = self._work_with_section(child, paragraphs, tag)
        return paragraphs

    def _extract_citations(self, tree, paragraphs):
        citations = {'sentence_ids': list(), 'papers': dict()}
        # find all labels of citations in text
        rid = tree.findall(".//xref[@ref-type='bibr']") + tree.findall(".//xref[@ref-type='other']") 
        rid_labels = {ref.attrib.get("rid", ''): "".join(ref.itertext()) for ref in rid if
                      ref.text is not None}
        for ref in tree.findall(".//ref"):
            rid = ref.attrib.get("id", '')
            # write an actual label (from the text)
            label = rid_labels.get(rid, getattr(ref.find('.//label'), 'text', ""))
            cit = {
                "text": [],
                "full_title": getattr(ref.find('.//article-title'), 'text', ""),
                "abstract": "",
                "journal": getattr(ref.find('.//source'), 'text', ""),
                "pmid": "",
                "pmc": "",
                "doi": "",
                "publisher_id": "",
                "author_list": [],
                "affiliation_list": [],
                "publication_year": getattr(ref.find('.//year'), 'text', ""),
                "publication_date": "",
                "subjects": []
            }
            pub_id = ref.find('.//pub-id')
            if pub_id is not None:
                cit[pub_id.attrib['pub-id-type']] = getattr(pub_id, 'text', "")
            citations['papers'][label] = cit
        # 4 in a row, wow
        for i, section in enumerate(paragraphs):
            for j, part in enumerate(section['section_text']):
                for k, sentence in enumerate(part):
                    for label in citations['papers']:
                        if label in sentence:
                            citations['sentence_ids'].append([i, j, k, label])
        return citations

    def __call__(self, data, data_format="pmc", parse_citations=False):
        """
        If parse_citations, then parse the paper and add citations section. Else - don't parse citations section.
        """
        paragraphs = OrderedDefaultdict(list)
        if data_format == "pmc" or data_format == "pubmed":
            parsed = self._parse_pmc_xml(data, paragraphs, parse_citations)
        else:
            raise NotImplementedError(f"The data format `{data_format}` is not yet supported by the parser.")
        return parsed

    def _parse_pmc_xml(self, data, paragraphs, parse_citations=True):
        '''
        # TODO: maybe we should deal with
        :param data: path to xml paper
        :param paragraphs:  OrderedDefaultdict
        :param parse_citations:
        :return:
        '''
        output = defaultdict(str)
        # load data from string
        tree = ET.fromstring(data)
        # extract abstract & text
        abstract = tree.find(".//abstract")
        if abstract:
            for ch in abstract.getchildren():
                if len(ch.text) > 50:
                    abstract = ch.text
                    break

        output['abstract'] = abstract

        sections = filter(None, [tree.find(".//body"), tree.find(".//floats-group")])
        for section in sections:
            paragraphs = self._work_with_section(section, paragraphs)
        # cast into appropriate format
        print(paragraphs)
        paragraphs = [
            {
                'section_title': key,
                'section_text': [
                    nltk.sent_tokenize(p) for p in paragraphs[key]
                ]
            }
            for key in paragraphs
        ]
        # del paragraphs['abstract']
        # print(paragraphs)
        output['text'] = paragraphs
        # : 'full_title',
        search_result = tree.find(".//" + 'article-title')
        if search_result is not None:
            title = "".join(search_result.itertext()).strip()
            output["full_title"] = " ".join(title.split()) # replace multiple spaces with one
        # extract full_title, journal, pmid, pmc, doi, publisher_id, publication_date, publication_year
        # keys are from xmls, values are from target json
        # maybe this argument should be included, 'publisher-name':'publisher_name'

        names_dict = {'journal-title': 'journal', 'pmid': 'pmid',
                      'pmc': 'pmc', 'doi': 'doi', 'publisher-id': 'publisher_id', 'copyright-year': 'publication_year'}
        for key in names_dict:
            search_result = tree.find(".//" + key)
            if search_result:
                output[names_dict[key]] = search_result.text
            search_results = tree.findall(".//article-id")
            for result in search_results:
                if result.attrib['pub-id-type'] == key:
                    output[names_dict[key]] = result.text

        # author_list
        author_list = []
        for author in tree.findall(".//contrib-group/contrib"):
            # surname, name, aff
            try:
                surname = author.find('./name/surname').text
            except AttributeError:
                pass
            try:
                name = author.find('./name/given-names').text
            except AttributeError:
                pass
            try:
                author_field = [[surname, name, 'aff' + p.text] for p in author.findall("./xref/[@ref-type='aff']/sup")]
            except AttributeError:
                continue
            if len(author_field):
                author_list.append(author_field)
        output['author_list'] = author_list

        # affiliation_list
        affiliation_list = []
        for aff in tree.findall(".//aff"):
            affiliation_field = []
            # 0 - for example 'aff1', 1 - name of institure
            affiliation_field.append(aff.attrib.get('id', ""))
            try:
                affiliation_field.append(aff.find("./addr-line").text)
            except AttributeError:
                affiliation_field.append("")
            affiliation_list.append(affiliation_field)
        output['affiliation_list'] = affiliation_list

        # parsing of citations
        if parse_citations:
            output['citations'] = self._extract_citations(tree, paragraphs)

        return output
