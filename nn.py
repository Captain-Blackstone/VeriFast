from sentence_transformers import SentenceTransformer, util


def dedup(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def get_sections_and_sentence(parsed_paper, sentence_id):
    section_id, paragraph_id, sentence_id, cit_label = sentence_id
    sentence = parsed_paper['text'][section_id][paragraph_id][sentence_id]
    sections = parsed_paper['citations']['papers'][cit_label]['text']
    return sentence, sections, cit_label


class SemanticSearcher:
    def __init__(self, model='msmarco-distilbert-base-tas-b', device='cpu', top_k=1):
        self.top_k = top_k
        self.model = SentenceTransformer('msmarco-distilbert-base-tas-b', device=device)

    def dirty_call(self, sentence, json_article, final_result):
        query_embedding = self.model.encode(sentence)
        sections = json_article["text"]
        main_sections = dedup([section["section_title"].split("|")[0] for section in sections])
        for main_section in main_sections:
            pps = []
            for section in sections:
                section_title, section_paragraphs = section["section_title"], section["section_text"]
                if main_section in section_title:
                    section_paragraphs = [" ".join(sentence_list) for sentence_list in section_paragraphs]
                    pps.extend(section_paragraphs)
            passage_embeddings = self.model.encode(pps)
            result = util.semantic_search(query_embedding, passage_embeddings, top_k=self.top_k,
                                          score_function=util.dot_score)
            result = pps[result[0][0]["corpus_id"]]
            final_result[main_section] = result
        return final_result


    def __call__(self, parsed_paper):
        # 1. add `selected_paragraphs` as a list of ids of citation sentences
        # 2. for each paper, for each section - run semantic search
        # 3. for each search result

        for sentence_id in parsed_paper['citations']['sentence_ids']:
            sentence, sections, cit_label = get_sections_and_sentence(parsed_paper, sentence_id)
            if 'selected_paragraphs' not in parsed_paper['citations']['papers'][cit_label]:
                parsed_paper['citations']['papers'][cit_label]['selected_paragraphs'] = [sentence_id]
            else:
                parsed_paper['citations']['papers'][cit_label]['selected_paragraphs'].append(sentence_id)
        for cit_label, paper in parsed_paper['citations']['papers'].values():
            sentences = [
                parsed_paper['text'][section_id][paragraph_id][sentence_id]
                for (section_id, paragraph_id, sentence_id, _) in paper['selected_paragraphs']
            ]
            sections = parsed_paper['citations']['papers'][cit_label]['text']
            query_embedding = self.model.encode(sentences, convert_to_tensor=True)
            selected_paragraphs = []
            for section in sections:
                paragraphs = ["".join(paragraph) for paragraph in section['section_text']]
                passage_embeddings = self.model.encode(paragraphs, convert_to_tensor=True)
                search_results = util.semantic_search(
                    query_embedding,
                    passage_embeddings,
                    top_k=self.top_k,
                    score_function=util.dot_score
                )
                selected_paragraphs.append(search_results)
            for i, element in enumerate(paper['selected_paragraphs']):
                element.append([(j, res) for j, res in enumerate(search_results)])
        return paper
