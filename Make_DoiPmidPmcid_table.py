from app import db#, text
from app.models import DoiPmcidPmid

def make_table_from_csv():
    for i, line in enumerate(text[1:]):
        if i % 100000 == 0:
            db.session.flush()
            print(i)
        try:
            pmid, pmcid, doi = line.strip().split(",")
        except:
            print(line)
            continue
        new_record = DoiPmcidPmid()
        if pmid:
            pmid = int(pmid)
            new_record.pmid = pmid
        if pmcid:
            new_record.pmcid = pmcid
        if doi != '""':
            new_record.doi = doi.replace("https://", "")
        db.session.add(new_record)
    db.session.commit()

def convert_dois_to_uniform_format():
    lastID = 0
    while True:
        things = db.session.query(DoiPmcidPmid).filter(DoiPmcidPmid.id > lastID).limit(1000).all()
        if not things or len(things) == 0: 
            break
        for article in things:
            lastID = article.id
            if article.doi is not None:
                article.doi = article.doi.split("doi.org/")[-1].replace('"', '')
            if lastID % 100000 == 0:
                db.session.flush()
                print(lastID, article.doi)
    
#    for i, article in enumerate(db.session.query(DoiPmcidPmid).all()):
#        article.doi = article.doi.split("doi.org/")[-1].replace('"', '')
#        if i % 100000 == 0:
#            db.session.flush()
#            print(i, article.doi)
    db.session.commit()


convert_dois_to_uniform_format()
