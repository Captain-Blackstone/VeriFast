"""articles table

Revision ID: 871b637c4dca
Revises: 
Create Date: 2021-09-24 10:05:35.878466

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '871b637c4dca'
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    op.create_table('article',
    sa.Column('id', sa.Integer(), nullable=False),
    sa.Column('doi', sa.String(length=32), nullable=True),
    sa.Column('title', sa.String(length=1000), nullable=True),
    sa.Column('filename', sa.String(length=50), nullable=True),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_index(op.f('ix_article_doi'), 'article', ['doi'], unique=True)
    op.create_index(op.f('ix_article_filename'), 'article', ['filename'], unique=True)
    op.create_index(op.f('ix_article_title'), 'article', ['title'], unique=True)
    # ### end Alembic commands ###


def downgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    op.drop_index(op.f('ix_article_title'), table_name='article')
    op.drop_index(op.f('ix_article_filename'), table_name='article')
    op.drop_index(op.f('ix_article_doi'), table_name='article')
    op.drop_table('article')
    # ### end Alembic commands ###
