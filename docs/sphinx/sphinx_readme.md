# Sphinx documentation for python-electric

This folder contains the **Sphinx-based HTML documentation setup** for the `python-electric` project.

The goal of this README is to serve as a **long-term reminder / checklist** explaining how the documentation is generated, updated, and maintained.

---

## Project context

- Project uses a **src-layout**:
  ```
  src/python_electric/
  ```
- A general `docs/` folder already exists and may contain:
  - drawings
  - examples
  - pdf
  - tryout
- Sphinx lives **isolated** in:
  ```
  docs/sphinx/
  ```

This avoids conflicts with other documentation artefacts.

---

## Folder structure

```
docs/sphinx/
├── Makefile          # Linux / macOS build helper
├── make.bat          # Windows build helper
├── source/
│   ├── conf.py       # Sphinx configuration
│   ├── index.rst     # Main documentation entry point
│   └── api/          # Auto-generated API documentation
└── build/            # Generated HTML output (NOT committed)
```

---

## 1. One-time setup

### 1.1 Install dependencies (in virtual environment)

```bash
python -m pip install sphinx sphinx-rtd-theme sphinxcontrib-napoleon
```

---

### 1.2 Initialize Sphinx (already done)

Executed from `docs/`:

```bash
sphinx-quickstart sphinx
```

Key choices:
- Separate source and build directories: **Yes**
- Project name: `python-electric`
- Language: `en`

---

## 2. `conf.py` essentials (src-layout)

Because the project uses a **src-layout**, Sphinx must be told where to find the package.

In `docs/sphinx/source/conf.py`:

```python
import os
import sys

# Make src/ importable for autodoc
sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), "../../../src")
))
```

Enabled extensions:

```python
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
]
```

Theme:

```python
html_theme = "sphinx_rtd_theme"
```

Optional (only if imports fail during build):

```python
# autodoc_mock_imports = ["pint"]
```

---

## 3. Generate API documentation

Run **from project root**:

```bash
sphinx-apidoc -o docs/sphinx/source/api src/python_electric
```

Recommended flags when regenerating:

```bash
sphinx-apidoc -o docs/sphinx/source/api src/python_electric -f -e
```

- `-f` → overwrite existing `.rst`
- `-e` → one page per module

---

## 4. Link API docs to main index

In `docs/sphinx/source/index.rst`:

```rst
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   api/modules
```

---

## 5. Build HTML documentation

Always run from **project root**:

```bash
python -m sphinx -b html -T -E -v docs/sphinx/source docs/sphinx/build/html
```

Output:

```
docs/sphinx/build/html/index.html
```

Open this file in a browser.

---

## 6. Daily workflow (short version)

After changing code or docstrings:

```bash
python -m sphinx -b html docs/sphinx/source docs/sphinx/build/html
```

---

## 7. Git rules

### Commit ✅

```
docs/sphinx/
  Makefile
  make.bat
  source/
```

### Do NOT commit ❌

```
docs/sphinx/build/
```

`.gitignore`:

```gitignore
docs/sphinx/build/
```

---

## 8. Common pitfalls (engineering projects)

- Import errors usually come from **side effects at import time**
- Always activate the correct virtual environment
- If Sphinx fails: rebuild with `-T -E -v`
- Do not document internal helper modules unless intentional

---

## 9. Recommended next extensions

- Add narrative pages:
  - `assumptions.rst`
  - `standards.rst`
  - `engineering_disclaimer.rst`
- Limit public API exposure
- Publish via ReadTheDocs or GitHub Pages

---

**This README is intended as a permanent reference for maintaining the Sphinx documentation of `python-electric`.**

