multiple_value_cols = (
    "standard_smiles",
    "canonical_smiles",
    "pchembl_value",
    "assay_chembl_id",
    "assay_description",
    "activity_id",
    "assay_type",
    "standard_type",
    "confidence_score",
    "standard_relation",
    "target_organism",
    "molecule_chembl_id",
    "document_chembl_id",
    "assay_tissue",
    "assay_cell_type",
    "relationship_type",
    "indication_class",
    "max_phase",
    "oral",
    "prodrug",
    "withdrawn_flag",
)

DEFAULT_ASSAY_MATCH_FIELDS = [
    "assay_type",
    "assay_organism",
    "assay_category",
    "assay_tax_id",
    "assay_strain",
    "assay_tissue",
    "assay_cell_type",
    "assay_subcellular_fraction",
    "bao_format",
]

DEFAULT_ASSAY_CHEMBL_ID = "assay_chembl_id"
DEFAULT_MOLECULE_CHEMBL_ID = "molecule_chembl_id"
DATA_DROPPING_COMMENT = "data_dropping_comment"
DATA_PROCESSING_COMMENT = "data_processing_comment"
