from app import db, text
from app.models import DoiPmcidPmid

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
