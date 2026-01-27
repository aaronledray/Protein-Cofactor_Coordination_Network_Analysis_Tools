# SSCNA Web App

A simple web interface for Single Structure Cofactor Network Analysis.

## Quick Start

```bash
# Install dependencies
pip install -r webapp/requirements.txt

# Run the app
streamlit run webapp/app.py
```

Then open http://localhost:8501 in your browser.

## Features

- Upload PDB or mmCIF structure files
- Specify cofactor residue names
- Adjust distance cutoffs
- View coordination sphere summaries
- Download results as CSV

## Deployment

### Streamlit Community Cloud (Free)

1. Push to GitHub
2. Go to [share.streamlit.io](https://share.streamlit.io)
3. Connect your repo and select `webapp/app.py`

### Hugging Face Spaces (Free)

1. Create a new Space with Streamlit SDK
2. Copy `webapp/app.py` and dependencies
3. Add `modules/` directory

### Docker

```dockerfile
FROM python:3.10-slim
WORKDIR /app
COPY . .
RUN pip install -r webapp/requirements.txt
EXPOSE 8501
CMD ["streamlit", "run", "webapp/app.py", "--server.address", "0.0.0.0"]
```
