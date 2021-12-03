import xml.etree.ElementTree as ET
from collections import defaultdict
import collections
from itertools import chain
import numpy as np
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

                tt = child.text
                if tt is None:
                    if len(list(child)) > 0:
                        tt = child[0].text
                if tt is not None:
                    if not tag:
                        tag = tt
                    else:
                        tag += "|" + tt

            if child.tag == "p":
                specital_paragraph = False
                paragraph_children = list(child)
                for p_child in paragraph_children:
                    if p_child.tag in ("fig", "table-wrap"):
                        paragraphs = self._work_with_section(p_child, paragraphs, tag)
                        specital_paragraph = True
                if not specital_paragraph:
                    paragraph = "".join(child.itertext()).strip()
                else:
                    parts = ([child.text] +
                             list(chain(*([c.text, c.tail] for c in list(child) if
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

    def _parse_xml(self, data, paragraphs, data_format, parse_citations=True):
        """
        :param data: path to xml paper
        :param paragraphs:  OrderedDefaultdict
        :param parse_citations:
        :return:
        """
        tree = ET.fromstring(data)
        if data_format in ["pmc", "pubmed"]:
            args = [tree, "pmc"]
        elif data_format == "elsevier":
            namespaces = {"bk": "http://www.elsevier.com/xml/bk/dtd",
                          "cals": "http://www.elsevier.com/xml/common/cals/dtd",
                          "ce": "http://www.elsevier.com/xml/common/dtd",
                          "ja": "http://www.elsevier.com/xml/ja/dtd",
                          "mml": "http://www.w3.org/1998/Math/MathML",
                          "sa": "http://www.elsevier.com/xml/common/struct-aff/dtd",
                          "sb": "http://www.elsevier.com/xml/common/struct-bib/dtd",
                          "tb": "http://www.elsevier.com/xml/common/table/dtd",
                          "xlink": "http://www.w3.org/1999/xlink",
                          "xocs": "http://www.elsevier.com/xml/xocs/dtd",
                          "dc": "http://purl.org/dc/elements/1.1/",
                          "dcterms": "http://purl.org/dc/terms/",
                          "prism": "http://prismstandard.org/namespaces/basic/2.0/",
                          "xsi": "http://www.w3.org/2001/XMLSchema-instance"}
            args = [tree, "elsevier", namespaces]

        output = defaultdict(str)
        output['abstract'] = self._find_abstract(*args)
        output['text'] = self._find_text(*args)
        output["full_title"] = self._find_full_title(*args)
        output['author_list'] = self._find_author_list(*args)
        output['affiliation_list'] = self._find_affiliation_list(*args)
        output = self._find_small_things(output, *args)
        if parse_citations:
            output['citations'] = self._find_citations(tree, output['text'], *args[1:])
        return output

    def _find_abstract(self, tree, publisher, namespaces=None):
        if publisher == "elsevier":
            abstract = tree.find(".//ce:abstract", namespaces)
            if abstract:
                abstract = "".join(abstract.itertext()).strip()
        elif publisher == "pmc":
            abstract = tree.find(".//abstract")
            if abstract:
                for ch in list(abstract):
                    if len(ch.text) > 50:
                        abstract = ch.text
                        break
        return abstract

    def _find_text(self, tree, publisher, namespaces=None):
        if publisher == "pmc":
            sections = filter(None, [tree.find(".//body"), tree.find(".//floats-group")])
            for section in sections:
                paragraphs = self._work_with_section(section, paragraphs)

            # cast into appropriate format
            paragraphs = [
                            {
                                'section_title': key,
                                'section_text': [nltk.sent_tokenize(p) for p in paragraphs[key]]
                            }
                            for key in paragraphs
                        ]
        elif publisher == "elsevier":
            raw_text = tree.find(".//xocs:rawtext", namespaces)
            structured_text = [el.strip() for el in raw_text.text.split("   ")]

            possible_megasections = ["Introduction", "Materials and Methods", "Results", "Discussion", "Conclusions",
                                     "Acknowledgements", "References"]
            sections_and_megasections = list(filter(lambda el: len(el.split()) < 20, structured_text))
            is_megasection = lambda title: any([title.lower() in el.lower() for el in possible_megasections])
            megasection_mask = [is_megasection(el) for el in sections_and_megasections]
            megasections_indices = [i for i in range(len(megasection_mask)) if megasection_mask[i]]
            sections_and_megasections = sections_and_megasections[min(megasections_indices):]
            megasection_mask = [is_megasection(el) for el in sections_and_megasections]
            megasections_indices = [i for i in range(len(megasection_mask)) if megasection_mask[i]]

            def find_corresponding_megasection(sct):
                deltas = np.array(megasections_indices) - sections_and_megasections.index(sct)
                megasection = sections_and_megasections[np.array(megasections_indices)[(deltas <= 0)][-1]]
                if sct != megasection:
                    return megasection + "|"
                else:
                    return ""

            sections_and_megasections_dict = {}
            for s in sections_and_megasections:
                sections_and_megasections_dict[s] = find_corresponding_megasection(s) + s
            paragraphs = []
            next_section = None
            for el in structured_text:
                if el in sections_and_megasections:
                    if next_section is not None and next_section["section_text"]:
                        paragraphs.append(next_section)
                    next_section = dict()
                    next_section["section_title"] = sections_and_megasections_dict[el]
                    next_section["section_text"] = []
                else:
                    if next_section is not None:
                        next_section["section_text"].append(nltk.sent_tokenize(el))
        return paragraphs

    def _find_full_title(self, tree, publisher, namespaces=None):
        if publisher == "pmc":
            search_result = tree.find(".//" + 'article-title')
            if search_result is not None:
                title = "".join(search_result.itertext()).strip()
                title = " ".join(title.split())  # replace multiple spaces with one
        elif publisher == "elsevier":
            search_result = tree.find(".//dc:title", namespaces)
            if search_result is not None:
                title = "".join(search_result.itertext()).strip()
                title = " ".join(title.split())  # replace multiple spaces with one
        return title

    def _find_author_list(self, tree, publisher, namespaces=None):
        author_list = []
        if publisher == "pmc":
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
                    author_field = [[surname, name, 'aff' + p.text] for p in
                                    author.findall("./xref/[@ref-type='aff']/sup")]
                except AttributeError:
                    continue
                if len(author_field):
                    author_list.append(author_field)
        elif publisher == "elsevier":
            for author in tree.findall(".//dc:creator", namespaces):
                author_list.append([author.text.split(", ") + [None]])
        return author_list

    def _find_affiliation_list(self, tree, publisher, namespaces=None):
        if publisher == "pmc":
            affiliation_list = []
            for aff in tree.findall(".//aff"):
                affiliation_field = [aff.attrib.get('id', "")]
                # 0 - for example 'aff1', 1 - name of institure
                try:
                    affiliation_field.append(aff.find("./addr-line").text)
                except AttributeError:
                    affiliation_field.append("")
                affiliation_list.append(affiliation_field)
        elif publisher == "elsevier":
            affiliation_list = []

        return affiliation_list

    def _find_small_things(self, output, tree, publisher, namespaces=None):
        if publisher == "pmc":
            # keys are from target json, values are from xmls
            names_dict = {'journal': 'journal-title', 'pmid': 'pmid',
                          'pmc': 'pmc', 'doi': 'doi', 'publisher_id': 'publisher-id',
                          'publication_year': 'copyright-year'}
            for key, val in names_dict.items():
                search_result = tree.find(".//" + val)
                if search_result:
                    output[key] = search_result.text
                search_results = tree.findall(".//article-id")
                for result in search_results:
                    if result.attrib['pub-id-type'] == val:
                        output[key] = result.text
        elif publisher == "elsevier":
            # keys are from target json, values are from xmls
            names_dict = {
                'journal': '-', 'pmid': '-', 'pmc': '-',
                'doi': 'prism:doi',
                'publisher_id': 'prism:publisher',
                'publication_year': '-'}
            for key, val in names_dict.items():
                search_result = tree.find(".//" + val, namespaces)
                if search_result is not None:
                    output[key] = "".join(search_result.itertext())
                else:
                    output[key] = None
        return output

    def _find_citations(self, tree, paragraphs, publisher, namespaces=None):
        citations = {'sentence_ids': list(), 'papers': dict()}
        if publisher == "pmc":
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
        elif publisher == "elsevier":
            references = tree.findall(".//ce:bib-reference", namespaces)
            for ref in references:
                label = ref.find("./ce:label", namespaces).text
                cit = {
                    "text": [],
                    "full_title": getattr(ref.find('./sb:reference/sb:contribution//sb:maintitle', namespaces), 'text', ""),
                    "abstract": "",
                    "journal": getattr(ref.find('./sb:reference/sb:host//sb:maintitle', namespaces), 'text', ""),
                    "pmid": "",
                    "pmc": "",
                    "doi": "",
                    "publisher_id": "",
                    "author_list": [],
                    "affiliation_list": [],
                    "publication_year": getattr(ref.find('./sb:reference/sb:host//sb:date', namespaces), 'text', ""),
                    "publication_date": "",
                    "subjects": []
                }
                citations["papers"][label] = cit
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
        if data_format in ["pmc", "elsevier", "pubmed"]:
            parsed = self._parse_xml(data, paragraphs, data_format, parse_citations)
        else:
            raise NotImplementedError(f"The data format `{data_format}` is not yet supported by the parser.")
        return parsed


if __name__ == '__main__':
    ParserObj = Parser()
    with open("/home/dmitry/Downloads/Telegram Desktop/elsevier_written_1.xml", "r", encoding="ISO-8859-1") as fl:
        xml = fl.read()
    text = ParserObj(xml, "elsevier", parse_citations=True)
    print(text)
