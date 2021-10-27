import json
import re
import xml.etree.ElementTree as ET
from habanero import Crossref
from Bio import Entrez


class Searcher:
    def __init__(self,
                 doi_pmid_pmcid_db,
                 api_key,
                 pubmed_email,
                 elsevier_apikey_config,
                 pmc_db="pmc",
                 pubmed_db="pubmed",
                 doi_pmcid_pmid_class=None):
        """
        This is a Searcher module. It is used to search for papers in the publishers API.
        The logic of this module is to start searching in one db,
        try to find at least meta-info (and use it in search), check whether the info is right (doi),
        proceed to another db until we get the text or we are out of options.
        Args:
            - pubmed_email: str - pubmed email used to surf pubmed database
        """
        self.doi_pmid_pmcid_db = doi_pmid_pmcid_db
        self.DoiPmcidPmidClass = doi_pmcid_pmid_class
        self.elsevier_apikey_config = elsevier_apikey_config
        self.pubmed_email = pubmed_email
        self.api_key = api_key
        self.pubmed_db = pubmed_db
        self.pmc_db = pmc_db
        self.doi_regexp = re.compile('\b(10[.][0-9]{4,}(?:[.][0-9]+)*/(?:(?!["&\'<>])\S)+)\b')
        self.cr = Crossref()

    async def __call__(self,
                       text_input=None,
                       doi=None,
                       pmcid=None,
                       pmid=None):
        if doi:
            search_query = doi.strip()
        elif pmcid:
            search_query = pmcid.strip()
        elif pmid:
            search_query = pmid.strip()
        elif text_input:
            search_query = text_input.strip()
        else:
            raise Exception("No valid input to searcher is provided")
        # Search Crossref database for the paper metadata
        if doi:
            crossref_meta = self.fetch_meta(search_query)
        else:
            crossref_meta = self.crossref_search_meta(search_query)
        try:
            text_link = crossref_meta["link"][0]["URL"]
        except KeyError:
            print(f"Could not get link for paper {search_query} from crossref: publisher is ueban")
            text_link = crossref_meta["URL"]
        publisher = crossref_meta["publisher"]
        id_tuple = self.convert_id(crossref_meta["DOI"])
        try:
            if id_tuple[1] == "":
                raise KeyError
            if "Elsevier" in publisher:
                raise KeyError
            pmc_xml = await self.fetch_fulltext(pmcid=id_tuple[1])
            pmc_parsed = ET.fromstring(pmc_xml)
            for child in pmc_parsed[0]:
                if child.tag == "body":
                    return pmc_xml, "pmc"
            raise KeyError
        except KeyError:
            try:
                pubmed_xml = self.ncbi_fetch(id_input=id_tuple[0], database="pubmed")
                return pubmed_xml, "pubmed"
            except:
                paper_xml = self.fetch_fulltext(url=text_link, publisher=publisher)
                return paper_xml, publisher

    def ncbi_search_id(self, search_query, database="pubmed"):
        """
        Return list with top PMC IDs by query and download papers from Id list
        """
        # getting search results for the query
        handle = Entrez.esearch(
            db=database,
            term=search_query + "[title]",
            prefix="xlink",
            sort="title",
        )
        search_results = Entrez.read(handle)
        handle.close()
        id_list = search_results["IdList"]  # list with PMC IDs
        if not len(id_list):
            print(f"No papers were found in pubmed for query {search_query}")
        else:
            return id_list

    def convert_id(self, doi):
        """
        Retrieve pmid and pmcid from doi of an article
        :param doi: doi of an article
        :return:
        """
        doi = doi.split("doi.org/")[-1].replace('"', '')
        record = self.doi_pmid_pmcid_db.session.query(self.DoiPmcidPmidClass).\
            filter(self.DoiPmcidPmidClass.doi == doi).first()
        print(record.pmid, record.pmcid, record.doi)
        return record.pmid, record.pmcid, record.doi

    def ncbi_fetch(self, id_input, database="pubmed"):
        # take a paper xml
        handle = Entrez.efetch(db=database, id=id_input, rettype="full", retmode="xml")
        # save xml if needed
        xml_data = handle.read()
        return xml_data

    def crossref_search_meta(self, query=None):
        """
        Requests a search in crossref database by paper name. Returns a json object with meta from the best match.
        :param query: a text query with paper name
        :return:
        """

        meta_json_search_results = self.cr.works(query=query)
        for item in meta_json_search_results["message"]["items"]:
            if re.match(item["title"][0], query, re.IGNORECASE):
                return item
        print(f"No exact match found for query {query}, returning closest search result")
        return meta_json_search_results["message"]["items"][0]

    async def fetch_meta(self,
                         doi=None,
                         pmid=None,
                         pmid_list=None,
                         database="crossref",
                         session=None):
        """
        Fetches paper metadata in json format by doi from crossref database
        """
        if database == "crossref":
            meta_json = self.cr.works(ids=doi)
            return meta_json["message"]
        elif database == "pubmed":
            if pmid:
                url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={pmid}&retmode=json"
                meta = await session.request(method="GET",
                                             url=url)
                return meta.json()
            elif pmid_list:
                url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/" \
                      f"esummary.fcgi?db=pubmed&id={','.join(pmid_list)}&retmode=json"
                meta = await session.request(method="GET",
                                             url=url)
                return meta.json()

    async def fetch_fulltext(self,
                             url=None,
                             pmcid=None,
                             publisher=None,
                             session=None):
        if url:
            # Place your kostily here:
            # For PLoS
            headers = {}
            if "PLoS" in publisher:
                url = (url + "&type=manuscript").replace("dx.plos.org/", "journals.plos.org/plosone/article/file?id=")
            if "Elsevier" in publisher:
                headers = {"X-ELS-APIKey": self.elsevier_apikey_config}
            publisher_response = await session.request(method="GET",
                                                       url=url,
                                                       headers=headers)
            if publisher_response.status == 200:
                text = await publisher_response.text()
                return text
            else:
                print(f"Could not load paper from url: {url}\n"
                      f"{publisher_response.status}: {publisher_response.reason}\n{publisher_response.url}")

        elif pmcid:
            pmc_data = self.ncbi_fetch(pmcid, database="pmc")
            return pmc_data


if __name__ == '__main__':
    from app import db
    from app.models import DoiPmcidPmid
    import json
    with open("config.json", "r") as config:
        config_dict = json.load(config)
        pubmed_email_config = config_dict["pubmed_email"]
        elsevier_apikey_configuration = config_dict["elsevier_apikey"]

    Entrez.email = pubmed_email_config
    searcher = Searcher(doi_pmid_pmcid_db=db,
                        doi_pmcid_pmid_class=DoiPmcidPmid,
                        api_key="",
                        pubmed_email=pubmed_email_config,
                        elsevier_apikey_config=elsevier_apikey_configuration,
                        pmc_db="pmc",
                        pubmed_db="pubmed")

    async def getres(title):
        print("starting getres")
        r1, r2 = await searcher(title)
        return r1
    print("1")
    res =  getres("Identification of SOX2 as a novel glioma-associated antigen and potential target "
                 "for T cell-based immunotherapy")
    print(res)
