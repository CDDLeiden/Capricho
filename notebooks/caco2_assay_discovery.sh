#!/bin/bash
# Discover and classify permeability assays by transport direction
# Covers Caco-2 and MDCK-MDR1 cell lines

# This search is more robust than simply '%apical to basolateral%' which misses
# several naming patterns in ChEMBL assay descriptions:
#   - "apical side to basolateral side" (extra "side" after apical breaks the match)
#   - "A to B" / "B to A" shorthand notation
#   - "basolateral side to apical side" (same issue for B→A)
#   - bidirectional assays ("apical to basolateral or basolateral to apical")
#   - averaged assays ("(PA-B+PB-A)/2")
#   - directionless permeability assays with no direction keyword at all

echo "================================================================"
echo "Permeability Assay Discovery & Classification"
echo "Cell lines: Caco-2, MDCK-MDR1"
echo "Generated: $(date)"
echo "================================================================"
echo ""

echo "################################################################"
echo "#                        CACO-2                                #"
echo "################################################################"
echo ""

# --------------------------------------------------------------------------
# 1. A→B assays — ALL known description patterns
#    → ChEMBL-all-caco2-a2b-assays.txt (used for A→B comparability analysis)
# --------------------------------------------------------------------------
echo "================================================================"
echo "1. A→B (Apical to Basolateral) Assays — Absorption"
echo "   Patterns: 'apical to basolateral', 'apical to basal',"
echo "             'apical side to basolateral side', 'A to B'"
echo "================================================================"

capricho explore --query "
SELECT a.chembl_id AS assay_chembl_id,
       SUBSTR(a.description, 1, 200) AS description,
       COUNT(DISTINCT act.molregno) AS n_compounds
FROM assays a
JOIN activities act ON a.assay_id = act.assay_id
WHERE LOWER(a.description) LIKE '%caco%'
  AND act.standard_units IS NOT NULL
  AND (   LOWER(a.description) LIKE '%apical to basolateral%'
       OR LOWER(a.description) LIKE '%apical to basal %'
       OR LOWER(a.description) LIKE '%apical side to basolateral%'
       OR LOWER(a.description) LIKE '%apical side to basal %'
       OR LOWER(a.description) LIKE '%apical (ph%to basolateral%'
       OR LOWER(a.description) LIKE '%apical (ph%to basal %'
       OR LOWER(a.description) LIKE '%(a to b)%'
       OR LOWER(a.description) LIKE '%in a to b direction%'
       OR LOWER(a.description) LIKE '%a-to-b%'
       OR LOWER(a.description) LIKE '%a to b direction%'
  )
  AND NOT (   LOWER(a.description) LIKE '%basolateral to apical%'
           OR LOWER(a.description) LIKE '%basal to apical%'
           OR LOWER(a.description) LIKE '%basolateral side to apical%'
           OR LOWER(a.description) LIKE '%b to a%'
  )
GROUP BY a.chembl_id, a.description
HAVING COUNT(DISTINCT act.molregno) >= 5
ORDER BY n_compounds DESC
"

echo ""

# --------------------------------------------------------------------------
# 2. B→A assays — ALL known description patterns
#    → ChEMBL-all-caco2-b2a-assays.txt (used for B→A comparability analysis)
# --------------------------------------------------------------------------
echo "================================================================"
echo "2. B→A (Basolateral to Apical) Assays — Efflux"
echo "   Patterns: 'basolateral to apical', 'basal to apical',"
echo "             'basolateral side to apical side', 'B to A'"
echo "================================================================"

capricho explore --query "
SELECT a.chembl_id AS assay_chembl_id,
       SUBSTR(a.description, 1, 200) AS description,
       COUNT(DISTINCT act.molregno) AS n_compounds
FROM assays a
JOIN activities act ON a.assay_id = act.assay_id
WHERE LOWER(a.description) LIKE '%caco%'
  AND act.standard_units IS NOT NULL
  AND (   LOWER(a.description) LIKE '%basolateral to apical%'
       OR LOWER(a.description) LIKE '%basal to apical%'
       OR LOWER(a.description) LIKE '%basolateral side to apical%'
       OR LOWER(a.description) LIKE '%basal side to apical%'
       OR LOWER(a.description) LIKE '%basolateral (ph%to apical%'
       OR LOWER(a.description) LIKE '%(b to a)%'
       OR LOWER(a.description) LIKE '%in b to a direction%'
       OR LOWER(a.description) LIKE '%b-to-a%'
       OR LOWER(a.description) LIKE '%b to a direction%'
       OR LOWER(a.description) LIKE '%basolateral permeab%'
  )
  AND NOT (   LOWER(a.description) LIKE '%apical to basolateral%'
           OR LOWER(a.description) LIKE '%apical to basal%'
           OR LOWER(a.description) LIKE '%apical side to basolateral%'
           OR LOWER(a.description) LIKE '%a to b%'
  )
GROUP BY a.chembl_id, a.description
HAVING COUNT(DISTINCT act.molregno) >= 5
ORDER BY n_compounds DESC
"

echo ""

# --------------------------------------------------------------------------
# 3. Bidirectional assays — mention BOTH directions or report averages/ratios
#    → included in ChEMBL-all-caco2-assays.txt (used for unit/type heterogeneity)
# --------------------------------------------------------------------------
echo "================================================================"
echo "3. Bidirectional / Averaged Assays"
echo "   These mention both directions or report combined metrics"
echo "================================================================"

capricho explore --query "
SELECT a.chembl_id AS assay_chembl_id,
       SUBSTR(a.description, 1, 250) AS description,
       COUNT(DISTINCT act.molregno) AS n_compounds
FROM assays a
JOIN activities act ON a.assay_id = act.assay_id
WHERE LOWER(a.description) LIKE '%caco%'
  AND act.standard_units IS NOT NULL
  AND (
       (    (LOWER(a.description) LIKE '%apical to basolateral%' OR LOWER(a.description) LIKE '%apical side to basolateral%' OR LOWER(a.description) LIKE '%a to b%')
        AND (LOWER(a.description) LIKE '%basolateral to apical%' OR LOWER(a.description) LIKE '%basolateral side to apical%' OR LOWER(a.description) LIKE '%b to a%')
       )
       OR LOWER(a.description) LIKE '%pa-b%pb-a%'
       OR LOWER(a.description) LIKE '%pb-a%pa-b%'
       OR LOWER(a.description) LIKE '%bi-directional%'
       OR LOWER(a.description) LIKE '%bidirectional%'
       OR LOWER(a.description) LIKE '%efflux ratio%'
  )
GROUP BY a.chembl_id, a.description
HAVING COUNT(DISTINCT act.molregno) >= 5
ORDER BY n_compounds DESC
"

echo ""

# --------------------------------------------------------------------------
# 4. UNMATCHED — Caco-2 permeability assays with NO direction keyword
#    These are the ones most likely to be miscategorised or missed entirely.
#    They need manual inspection (e.g., checking the original paper).
#    → included in ChEMBL-all-caco2-assays.txt (used for unit/type heterogeneity)
# --------------------------------------------------------------------------
echo "================================================================"
echo "4. UNMATCHED — Permeability assays WITHOUT direction keywords"
echo "   These require manual inspection of the source literature"
echo "================================================================"

capricho explore --query "
SELECT a.chembl_id AS assay_chembl_id,
       SUBSTR(a.description, 1, 200) AS description,
       COUNT(DISTINCT act.molregno) AS n_compounds
FROM assays a
JOIN activities act ON a.assay_id = act.assay_id
WHERE LOWER(a.description) LIKE '%caco%'
  AND act.standard_units IS NOT NULL
  AND (   LOWER(a.description) LIKE '%permeab%'
       OR LOWER(a.description) LIKE '%papp%'
  )
  AND NOT (
          LOWER(a.description) LIKE '%apical to basolateral%'
       OR LOWER(a.description) LIKE '%apical to basal%'
       OR LOWER(a.description) LIKE '%apical side to basolateral%'
       OR LOWER(a.description) LIKE '%apical side to basal%'
       OR LOWER(a.description) LIKE '%basolateral to apical%'
       OR LOWER(a.description) LIKE '%basal to apical%'
       OR LOWER(a.description) LIKE '%basolateral side to apical%'
       OR LOWER(a.description) LIKE '%basal side to apical%'
       OR LOWER(a.description) LIKE '%(a to b)%'
       OR LOWER(a.description) LIKE '%(b to a)%'
       OR LOWER(a.description) LIKE '%a to b direction%'
       OR LOWER(a.description) LIKE '%b to a direction%'
       OR LOWER(a.description) LIKE '%apical (ph%to basolateral%'
       OR LOWER(a.description) LIKE '%basolateral (ph%to apical%'
  )
GROUP BY a.chembl_id, a.description
HAVING COUNT(DISTINCT act.molregno) >= 5
ORDER BY n_compounds DESC
"

echo ""

# --------------------------------------------------------------------------
# 5. Summary counts (not used directly in notebook)
# --------------------------------------------------------------------------
echo "================================================================"
echo "5. Summary: Assay counts by classification"
echo "================================================================"

capricho explore --query "
SELECT
  CASE
    WHEN (    (LOWER(a.description) LIKE '%apical to basolateral%' OR LOWER(a.description) LIKE '%apical to basal %' OR LOWER(a.description) LIKE '%apical side to basolateral%' OR LOWER(a.description) LIKE '%apical side to basal %')
          AND (LOWER(a.description) LIKE '%basolateral to apical%' OR LOWER(a.description) LIKE '%basal to apical%' OR LOWER(a.description) LIKE '%basolateral side to apical%')
    ) THEN 'bidirectional'
    WHEN LOWER(a.description) LIKE '%apical to basolateral%'
      OR LOWER(a.description) LIKE '%apical to basal %'
      OR LOWER(a.description) LIKE '%apical side to basolateral%'
      OR LOWER(a.description) LIKE '%apical side to basal %'
      THEN 'A_to_B'
    WHEN LOWER(a.description) LIKE '%basolateral to apical%'
      OR LOWER(a.description) LIKE '%basal to apical%'
      OR LOWER(a.description) LIKE '%basolateral side to apical%'
      OR LOWER(a.description) LIKE '%basal side to apical%'
      THEN 'B_to_A'
    ELSE 'no_direction'
  END AS direction,
  COUNT(DISTINCT a.chembl_id) AS n_assays,
  COUNT(DISTINCT act.molregno) AS n_compounds,
  COUNT(*) AS n_measurements
FROM assays a
JOIN activities act ON a.assay_id = act.assay_id
WHERE LOWER(a.description) LIKE '%caco%'
  AND act.standard_units IS NOT NULL
  AND (LOWER(a.description) LIKE '%permeab%' OR LOWER(a.description) LIKE '%papp%')
GROUP BY direction
ORDER BY n_measurements DESC
"

echo ""
echo "################################################################"
echo "#                       MDCK-MDR1                              #"
echo "################################################################"
echo ""

# --------------------------------------------------------------------------
# 6. MDCK-MDR1 A→B assays
#    → ChEMBL-all-mdck-mdr1-a2b-assays.txt (used for A→B comparability analysis)
# --------------------------------------------------------------------------
echo "================================================================"
echo "6. MDCK-MDR1 A→B (Apical to Basolateral) Assays"
echo "   Cell filter: description LIKE '%mdck%mdr1%' OR '%mdr1%mdck%'"
echo "   P-gp inhibitor assays INCLUDED: P-gp affects efflux (B→A),"
echo "   not absorption (A→B), so inhibitor presence is irrelevant"
echo "================================================================"

capricho explore --query "
SELECT a.chembl_id AS assay_chembl_id,
       SUBSTR(a.description, 1, 200) AS description,
       COUNT(DISTINCT act.molregno) AS n_compounds
FROM assays a
JOIN activities act ON a.assay_id = act.assay_id
WHERE (LOWER(a.description) LIKE '%mdck%mdr1%' OR LOWER(a.description) LIKE '%mdr1%mdck%')
  AND act.standard_units IS NOT NULL
  AND act.standard_units != '%'
  AND (act.standard_type IN ('Papp', 'permeability', 'logPapp', 'Pc', 'Pm')
       OR LOWER(act.standard_type) LIKE '%papp%')
  AND (   LOWER(a.description) LIKE '%apical to basolateral%'
       OR LOWER(a.description) LIKE '%apical to basal%'
       OR LOWER(a.description) LIKE '%apical side to basolateral%'
       OR LOWER(a.description) LIKE '%apical side to basal %'
       OR LOWER(a.description) LIKE '%apical (ph%to basolateral%'
       OR LOWER(a.description) LIKE '%apical (ph%to basal %'
       OR LOWER(a.description) LIKE '%(a to b)%'
       OR LOWER(a.description) LIKE '%(a)-to-basolateral%'
       OR LOWER(a.description) LIKE '%in a to b direction%'
       OR LOWER(a.description) LIKE '%a-to-b%'
       OR LOWER(a.description) LIKE '%a to b direction%'
  )
  AND NOT (   LOWER(a.description) LIKE '%basolateral to apical%'
           OR LOWER(a.description) LIKE '%basal to apical%'
           OR LOWER(a.description) LIKE '%basolateral side to apical%'
           OR LOWER(a.description) LIKE '%b to a%'
  )
GROUP BY a.chembl_id, a.description
HAVING COUNT(DISTINCT act.molregno) >= 5
ORDER BY n_compounds DESC
"

echo ""

# --------------------------------------------------------------------------
# 7. MDCK-MDR1 B→A assays
#    → ChEMBL-all-mdck-mdr1-b2a-assays.txt (used for B→A comparability analysis)
# --------------------------------------------------------------------------
echo "================================================================"
echo "7. MDCK-MDR1 B→A (Basolateral to Apical) Assays"
echo "   Excludes: assays performed in presence of P-gp inhibitors"
echo "================================================================"

capricho explore --query "
SELECT a.chembl_id AS assay_chembl_id,
       SUBSTR(a.description, 1, 200) AS description,
       COUNT(DISTINCT act.molregno) AS n_compounds
FROM assays a
JOIN activities act ON a.assay_id = act.assay_id
WHERE (LOWER(a.description) LIKE '%mdck%mdr1%' OR LOWER(a.description) LIKE '%mdr1%mdck%')
  AND act.standard_units IS NOT NULL
  AND act.standard_units != '%'
  AND (act.standard_type IN ('Papp', 'permeability', 'logPapp', 'Pc', 'Pm')
       OR LOWER(act.standard_type) LIKE '%papp%')
  AND (   LOWER(a.description) LIKE '%basolateral to apical%'
       OR LOWER(a.description) LIKE '%basal to apical%'
       OR LOWER(a.description) LIKE '%basolateral side to apical%'
       OR LOWER(a.description) LIKE '%basal side to apical%'
       OR LOWER(a.description) LIKE '%basolateral (ph%to apical%'
       OR LOWER(a.description) LIKE '%(b to a)%'
       OR LOWER(a.description) LIKE '%(b)-to-apical%'
       OR LOWER(a.description) LIKE '%in b to a direction%'
       OR LOWER(a.description) LIKE '%b-to-a%'
       OR LOWER(a.description) LIKE '%b to a direction%'
  )
  AND NOT (   LOWER(a.description) LIKE '%apical to basolateral%'
           OR LOWER(a.description) LIKE '%apical to basal%'
           OR LOWER(a.description) LIKE '%apical side to basolateral%'
           OR LOWER(a.description) LIKE '%a to b%'
  )
  AND NOT (   LOWER(a.description) LIKE '%inhibitor%'
           OR LOWER(a.description) LIKE '%in presence%'
           OR LOWER(a.description) LIKE '%in the presence%'
           OR LOWER(a.description) LIKE '%presence of%'
           OR LOWER(a.description) LIKE '%lsn335984%'
           OR LOWER(a.description) LIKE '%lsn 335984%'
           OR LOWER(a.description) LIKE '%gf120918%'
           OR LOWER(a.description) LIKE '%gf120919%'
           OR LOWER(a.description) LIKE '%elacridar%'
           OR LOWER(a.description) LIKE '%verapamil%'
           OR LOWER(a.description) LIKE '%cyclosporin%'
           OR LOWER(a.description) LIKE '%valspodar%'
           OR LOWER(a.description) LIKE '%zosuquidar%'
  )
GROUP BY a.chembl_id, a.description
HAVING COUNT(DISTINCT act.molregno) >= 5
ORDER BY n_compounds DESC
"

echo ""

# --------------------------------------------------------------------------
# 8. MDCK-MDR1 directionless permeability assays
#    → included in ChEMBL-all-mdck-mdr1-assays.txt (used for unit/type heterogeneity)
# --------------------------------------------------------------------------
echo "================================================================"
echo "8. MDCK-MDR1 — Permeability assays WITHOUT direction keywords"
echo "================================================================"

capricho explore --query "
SELECT a.chembl_id AS assay_chembl_id,
       SUBSTR(a.description, 1, 200) AS description,
       COUNT(DISTINCT act.molregno) AS n_compounds
FROM assays a
JOIN activities act ON a.assay_id = act.assay_id
WHERE (LOWER(a.description) LIKE '%mdck%mdr1%' OR LOWER(a.description) LIKE '%mdr1%mdck%')
  AND act.standard_units IS NOT NULL
  AND act.standard_units != '%'
  AND (act.standard_type IN ('Papp', 'permeability', 'logPapp', 'Pc', 'Pm')
       OR LOWER(act.standard_type) LIKE '%papp%')
  AND NOT (
          LOWER(a.description) LIKE '%apical to basolateral%'
       OR LOWER(a.description) LIKE '%apical to basal%'
       OR LOWER(a.description) LIKE '%apical side to basolateral%'
       OR LOWER(a.description) LIKE '%apical side to basal%'
       OR LOWER(a.description) LIKE '%basolateral to apical%'
       OR LOWER(a.description) LIKE '%basal to apical%'
       OR LOWER(a.description) LIKE '%basolateral side to apical%'
       OR LOWER(a.description) LIKE '%basal side to apical%'
       OR LOWER(a.description) LIKE '%(a to b)%'
       OR LOWER(a.description) LIKE '%(b to a)%'
       OR LOWER(a.description) LIKE '%a to b direction%'
       OR LOWER(a.description) LIKE '%b to a direction%'
       OR LOWER(a.description) LIKE '%apical (ph%to basolateral%'
       OR LOWER(a.description) LIKE '%basolateral (ph%to apical%'
  )
GROUP BY a.chembl_id, a.description
HAVING COUNT(DISTINCT act.molregno) >= 5
ORDER BY n_compounds DESC
"

echo ""
