#!/usr/bin/env bash
# ==============================================================================
# One-Click Automated Setup & Pipeline Execution for AlphaGenome PRNP Analysis
# ==============================================================================
set -e

echo "=========================================================================="
echo "  AlphaGenome (Google DeepMind, Nature 2026) PRNP Analysis Pipeline"
echo "=========================================================================="

# 1. Check Python version
if ! command -v python3 &>/dev/null; then
    echo "❌ Error: python3 is not installed. Please install Python 3.10+."
    exit 1
fi

PYTHON_VERSION=$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')
echo "✓ Detected Python version: $PYTHON_VERSION"

# 2. Virtual Environment Setup
VENV_DIR=".venv"
if [ ! -d "$VENV_DIR" ]; then
    echo "Creating virtual environment at $VENV_DIR..."
    python3 -m venv "$VENV_DIR"
fi

echo "Activating virtual environment..."
source "$VENV_DIR/bin/activate"

# 3. Install Requirements
echo "Installing dependencies from requirements.txt..."
pip install --upgrade pip --quiet
pip install -r requirements.txt --quiet
echo "✓ All dependencies installed successfully."

# 4. API Key Configuration
if [ -z "$ALPHAGENOME_API_KEY" ]; then
    if [ -f ".env" ]; then
        echo "Loading ALPHAGENOME_API_KEY from local .env file..."
        export $(grep -v '^#' .env | xargs)
    else
        echo "Note: ALPHAGENOME_API_KEY environment variable is not set."
        echo "To run live inference, set ALPHAGENOME_API_KEY in your shell or .env file."
    fi
fi

# 5. Execute Analysis & Figure Generation
echo ""
echo "Running Master Analysis Pipeline..."
python run_all.py

echo ""
echo "=========================================================================="
echo "  🎉 Pipeline Finished Successfully!"
echo "  - Clean Figures:      figures/clean/"
echo "  - Annotated Figures:  figures/annotated/"
echo "  - Reference Pages:    figures/paper_reference/"
echo "  - Slide Deck:         docs/slides_with_narration.md"
echo "=========================================================================="
