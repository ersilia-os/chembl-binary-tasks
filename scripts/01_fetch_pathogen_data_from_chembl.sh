SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

# Function to display usage
usage() {
    echo "Usage: $0 --pathogen_code <code> --output_dir <directory>"
    echo "  --pathogen_code   Specify the pathogen code (e.g., abaumannii)"
    echo "  --output_dir      Specify the output directory (absolute path)"
    exit 1
}

# Check for arguments
if [ $# -lt 4 ]; then
    echo "Error: Insufficient arguments."
    usage
fi

# Initialize variables
PATHOGEN_CODE=""
OUTPUT_DIR=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --pathogen_code)
            PATHOGEN_CODE=$2
            shift 2
            ;;
        --output_dir)
            OUTPUT_DIR=$2
            shift 2
            ;;
        *)
            echo "Error: Unknown argument $1"
            usage
            ;;
    esac
done

# Validate arguments
if [ -z "$PATHOGEN_CODE" ]; then
    echo "Error: --pathogen_code is required."
    usage
fi

if [ -z "$OUTPUT_DIR" ]; then
    echo "Error: --output_dir is required."
    usage
fi

echo "Fetching data for pathogen $PATHOGEN_CODE and storing in $OUTPUT_DIR"

python $SCRIPT_DIR/../src/pathogen_getter.py --pathogen_code $PATHOGEN_CODE --output_dir $OUTPUT_DIR 
python $SCRIPT_DIR/../src/clean_fetched_pathogen_data.py --pathogen_code $PATHOGEN_CODE --output_dir $OUTPUT_DIR
python $SCRIPT_DIR/../src/binarize_fetched_pathogen_data.py --pathogen_code $PATHOGEN_CODE --output_dir $OUTPUT_DIR
python $SCRIPT_DIR/../src/datasets_modelability.py --pathogen_code $PATHOGEN_CODE --output_dir $OUTPUT_DIR
python $SCRIPT_DIR/../src/select_tasks.py --pathogen_code $PATHOGEN_CODE --output_dir $OUTPUT_DIR
python $SCRIPT_DIR/../src/wrapup_tasks_and_clean_output_folder.py --pathogen_code $PATHOGEN_CODE --output_dir $OUTPUT_DIR