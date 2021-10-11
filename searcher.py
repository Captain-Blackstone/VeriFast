import json
import re
import xml.etree.ElementTree as ET
import urllib.request as libreq

import requests
from habanero import Crossref

from Bio import Entrez
from requests import get


class Searcher:
    def __init__(self, api_key, pubmed_email, elsevier_apikey_config, pmc_db="pmc", pubmed_db="pubmed"):
        """
        This is a Searcher module. It is used to search for papers in the publishers API.
        The logic of this module is to start searching in one db,
        try to find at least meta-info (and use it in search), check whether the info is right (doi), proceed to another db until we get the text or we are out of options.
        Args:
            - pubmed_email: str - pubmed email used to surf pubmed database
        """
        self.elsevier_apikey_config = elsevier_apikey_config
        self.pubmed_email = pubmed_email
        self.api_key = api_key
        self.pubmed_db = pubmed_db
        self.pmc_db = pmc_db
        self.doi_regexp = re.compile('\b(10[.][0-9]{4,}(?:[.][0-9]+)*/(?:(?!["&\'<>])\S)+)\b')
        self.cr = Crossref()

    def __call__(self, text_input=None, doi=None, pmcid=None, pmid=None):
        db = self.pubmed_db
        if doi:
            search_query = doi.strip()
        elif pmcid:
            search_query = pmcid.strip()
            db = self.pmc_db
        elif pmid:
            search_query = pmid.strip()
        elif text_input:
            search_query = text_input.strip()
        else:
            raise Exception("No valid input to searcher is provided")
        # Search Crossref database for the paper metadata
        print("s1")
        if doi:
            print("s1_doi")
            crossref_meta = self.fetch_meta(search_query)
        else:
            print("s1_no_doi")
            crossref_meta = self.crossref_search_meta(search_query)
        print("s2")
        try:
            text_link = crossref_meta["link"][0]["URL"]
        except KeyError:
            print(f"Could not get link for paper {search_query} from crossref: publisher is ueban")
            text_link = crossref_meta["URL"]
        publisher = crossref_meta["publisher"]
        id_tuple = self.convert_id(crossref_meta["DOI"])
        print("s3")
        try:
            print(id_tuple)
            if id_tuple[1] == "":
                raise KeyError
            if "Elsevier" in publisher:
                raise KeyError
            pmc_xml = self.fetch_fulltext(pmcid=id_tuple[1])
            pmc_parsed = ET.fromstring(pmc_xml)
            for child in pmc_parsed[0]:
                if child.tag == "body":
                    return pmc_xml, "pmc"
            raise KeyError
        except KeyError:
            print("s4e")
            try:
                pubmed_xml = self.ncbi_fetch(id_input=id_tuple[0], db="pubmed")
                return pubmed_xml, "pubmed"
            except:
                paper_xml = self.fetch_fulltext(url=text_link, publisher=publisher)
                return paper_xml, publisher


    def ncbi_search_id(self, search_query, db="pubmed", papers_count=20):
        """
        Return list with top PMC IDs by query and download papers from Id list
        """
        # getting search results for the query
        handle = Entrez.esearch(
            db=db,
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
        Convert a list of pmids or an individual pmid into a list of doi's
        :param id_list:
        :param convert_from:
        :param convert_to:
        :param db:
        :return:
        """
        record = db.session.query(DoiPmcidPmid).filter(DoiPmcidPmid.doi == doi).one()
        return record.pmid, record.pmcid, record.doi

    def ncbi_fetch(self, id_input, db="pubmed"):
        ### take a paper xml
        handle = Entrez.efetch(db=db, id=id_input, rettype="full", retmode="xml")
        # save xml if needed
        xml_data = handle.read()
        return xml_data

    def crossref_search_meta(self, query=None):
        """
        Requests a search in crossref database by paper name. Returns a json object with meta from the best match.
        :param query: a text query with paper name
        :return:
        """
        print("starting_cr_work")
        print(query)
        meta_json_search_results = self.cr.works(query=query)
        print("cr_worked")
        for item in meta_json_search_results["message"]["items"]:
            if re.match(item["title"][0], query, re.IGNORECASE):
                return item
        print(f"No exact match found for query {query}, returning closest search result")
        return meta_json_search_results["message"]["items"][0]

    def fetch_meta(self, doi=None, pmid=None, pmid_list=None, db="crossref"):
        """
        Fetches paper metadata in json format by doi from crossref database
        :param doi:
        :return:
        """
        if db == "crossref":
            meta_json = self.cr.works(ids=doi)
            return meta_json["message"]
        elif db == "pubmed":
            if pmid:
                meta_json = get(
                    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={}&retmode=json".format(
                        pmid)).json()
                return meta_json
            elif pmid_list:
                meta_json = get(
                    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={}&retmode=json".format(
                        ",".join(pmid_list))).json()
                return meta_json

    def fetch_fulltext(self, url=None, pmcid=None, publisher=None, email=None):
        if url:
            # Place your kostily here:
            # For PLoS
            if "PLoS" in publisher:
                url = (url + "&type=manuscript").replace("dx.plos.org/", "journals.plos.org/plosone"
                                                                                           "/article/file?id=")
            if "Elsevier" in publisher:
                headers = {"X-ELS-APIKey": self.elsevier_apikey_config}
            publisher_response = requests.get(url, headers=headers)
            if publisher_response.status_code == 200:
                return publisher_response.text
            else:
                print(f"Could not load paper from url: {url}\n"
                                f"{publisher_response.status_code}: {publisher_response.reason}\n{publisher_response.url}")

        elif pmcid:
            pmc_data = self.ncbi_fetch(pmcid, db="pmc")
            return pmc_data


if __name__ == '__main__':
    with open("config.json", "r") as config:
        config_dict = json.load(config)
        pubmed_email_config = config_dict["pubmed_email"]
        elsevier_apikey_config = config_dict["elsevier_apikey"]

    with open("PMID_PMCID_DOI.csv", "r") as fl:
        text = fl.readlines()

    Entrez.email = pubmed_email_config
    searcher = Searcher(api_key="", pubmed_email=pubmed_email_config, elsevier_apikey_config=elsevier_apikey_config,
                        pmid_pmcid_doi_db_text=text, pmc_db="pmc", pubmed_db="pubmed")
    print(searcher(
        "Genomic Profiles in Stage I Primary Non Small Cell Lung Cancer Using Comparative Genomic Hybridization Analysis of cDNA Microarrays"))
