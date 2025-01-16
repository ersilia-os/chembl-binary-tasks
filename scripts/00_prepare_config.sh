SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

echo "Preparing config file"
python $SCRIPT_DIR/../src/units_resolver_using_statistics.py

echo "Running LLM"
python $SCRIPT_DIR/../src/activity_binarizer_and_directions_with_llm.py

echo "Getting all molecules from ChEMBL for the molecule sampler"
python $SCRIPT_DIR/../src/chembl_chemical_space.py

echo "Cleaning activities"
python $SCRIPT_DIR/../src/clean_all_activities.py