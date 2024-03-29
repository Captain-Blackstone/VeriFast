"""DoiPmidPmcid table

Revision ID: 8bc4fddd9f95
Revises: 871b637c4dca
Create Date: 2021-09-27 15:17:36.523233

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '8bc4fddd9f95'
down_revision = '871b637c4dca'
branch_labels = None
depends_on = None


def upgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    op.create_table('doi_pmcid_pmid',
    sa.Column('id', sa.Integer(), nullable=False),
    sa.Column('doi', sa.String(length=64), nullable=True),
    sa.Column('pmid', sa.Integer(), nullable=True),
    sa.Column('pmcid', sa.Integer(), nullable=True),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_index(op.f('ix_doi_pmcid_pmid_doi'), 'doi_pmcid_pmid', ['doi'], unique=True)
    op.create_index(op.f('ix_doi_pmcid_pmid_pmcid'), 'doi_pmcid_pmid', ['pmcid'], unique=True)
    op.create_index(op.f('ix_doi_pmcid_pmid_pmid'), 'doi_pmcid_pmid', ['pmid'], unique=True)
    # ### end Alembic commands ###


def downgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    op.drop_index(op.f('ix_doi_pmcid_pmid_pmid'), table_name='doi_pmcid_pmid')
    op.drop_index(op.f('ix_doi_pmcid_pmid_pmcid'), table_name='doi_pmcid_pmid')
    op.drop_index(op.f('ix_doi_pmcid_pmid_doi'), table_name='doi_pmcid_pmid')
    op.drop_table('doi_pmcid_pmid')
    # ### end Alembic commands ###
