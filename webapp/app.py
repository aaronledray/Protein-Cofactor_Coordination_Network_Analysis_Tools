"""
SSCNA Web App - Single Structure Cofactor Network Analysis

A simple Streamlit interface for analyzing protein-cofactor coordination networks.

Run with: streamlit run webapp/app.py
"""

import os
import sys
import tempfile
import io

import streamlit as st
import pandas as pd

# Add project root to path
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from modules.io_utils import unpack_pdb_file
from modules.analysis import identify_coordination_network, generate_coordination_csv_with_moieties
from modules.moieties import bond_lookup
from modules.plotting import plot_interactive_modes_with_network
from modules.moieties import atom_type_colors

# Page config
st.set_page_config(
    page_title="SSCNA - Cofactor Network Analysis",
    layout="wide",
)

# Custom CSS for white background with readable text
st.markdown("""
<style>
    /* Main app background */
    .stApp {
        background-color: white;
        color: #1a1a1a;
    }

    /* Header/top bar */
    header[data-testid="stHeader"] {
        background-color: white;
    }

    /* Sidebar */
    .stSidebar, [data-testid="stSidebar"], [data-testid="stSidebarContent"] {
        background-color: #f8f9fa;
    }

    /* Text elements */
    .stMarkdown, .stText, p, span, label, .stMetric label, .stMetric [data-testid="stMetricValue"] {
        color: #1a1a1a !important;
    }
    h1, h2, h3, h4, h5, h6 {
        color: #1a1a1a !important;
    }
    .stTabs [data-baseweb="tab"] {
        color: #1a1a1a;
    }
    .stDataFrame {
        color: #1a1a1a;
    }

    /* File uploader */
    [data-testid="stFileUploader"], [data-testid="stFileUploaderDropzone"] {
        background-color: #f8f9fa !important;
        border-color: #dee2e6 !important;
    }
    [data-testid="stFileUploaderDropzone"] span, [data-testid="stFileUploaderDropzone"] small {
        color: #1a1a1a !important;
    }

    /* Text inputs */
    .stTextInput input, [data-testid="stTextInput"] input {
        background-color: white !important;
        color: #1a1a1a !important;
        border-color: #dee2e6 !important;
    }

    /* Slider */
    .stSlider label, .stSlider [data-testid="stTickBarMin"], .stSlider [data-testid="stTickBarMax"] {
        color: #1a1a1a !important;
    }

    /* Expander */
    .streamlit-expanderHeader {
        background-color: #f8f9fa !important;
        color: #1a1a1a !important;
    }

    /* Info box */
    .stAlert {
        background-color: #e7f3ff !important;
        color: #1a1a1a !important;
    }
</style>
""", unsafe_allow_html=True)

# Header
st.title("SSCNA: Single Structure Cofactor Network Analysis")
st.markdown("""
Analyze coordination networks around protein cofactors. Upload a structure file,
specify the cofactor, and visualize the primary and secondary coordination spheres.
""")

# Sidebar for inputs
with st.sidebar:
    st.header("Analysis Parameters")

    # File upload
    uploaded_file = st.file_uploader(
        "Upload Structure File",
        type=["pdb", "cif", "mmcif"],
        help="PDB or mmCIF format structure file"
    )

    # Cofactor input
    cofactor_input = st.text_input(
        "Cofactor Residue Name(s)",
        value="CU",
        help="Comma-separated residue names (e.g., 'CU' or 'ICS,CLF')"
    )

    # Optional second cofactor
    cofactor2_input = st.text_input(
        "Secondary Cofactor (optional)",
        value="",
        help="For combinatorial mode - residues that must be near primary cofactor"
    )

    # Distance cutoff
    distance_cutoff = st.slider(
        "Distance Cutoff (Å)",
        min_value=2.0,
        max_value=6.0,
        value=3.6,
        step=0.1,
        help="Maximum distance from cofactor for primary coordination sphere"
    )

    # Expand residues option
    expand_residues = st.checkbox(
        "Expand to Full Residues",
        value=False,
        help="Include all atoms of coordinating residues (not just seed atoms)"
    )

    # Combinatorial mode
    combinatorial = st.checkbox(
        "Combinatorial Mode",
        value=False,
        help="Filter cofactors by proximity to secondary cofactor"
    )

    if combinatorial:
        comb_cutoff = st.slider(
            "Combinatorial Cutoff (Å)",
            min_value=5.0,
            max_value=30.0,
            value=20.0,
            step=1.0,
        )
    else:
        comb_cutoff = 20.0

    # Run button
    run_analysis = st.button("Run Analysis", type="primary", use_container_width=True)

# Main content area
if uploaded_file is not None and run_analysis:
    # Parse inputs
    cofactor_list = [c.strip().upper() for c in cofactor_input.split(",") if c.strip()]
    cofactor2_list = [c.strip().upper() for c in cofactor2_input.split(",") if c.strip()] if cofactor2_input else None

    if not cofactor_list:
        st.error("Please enter at least one cofactor residue name.")
    else:
        with st.spinner("Analyzing structure..."):
            try:
                # Save uploaded file to temp location
                suffix = ".pdb" if uploaded_file.name.endswith(".pdb") else ".cif"
                with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
                    tmp.write(uploaded_file.getvalue())
                    tmp_path = tmp.name

                # Load structure
                structure, atoms_data = unpack_pdb_file(tmp_path)

                # Create temp output directory
                with tempfile.TemporaryDirectory() as output_dir:
                    # Run analysis
                    cofactor_sphere, pcs_atoms, scs_atoms = identify_coordination_network(
                        structure=structure,
                        cofactor_resname=cofactor_list,
                        distance_cutoff=distance_cutoff,
                        expand_residues=expand_residues,
                        combinatorial_mode=combinatorial,
                        combinatorial_cofactor_cutoff=comb_cutoff,
                        cofactor_resname2=cofactor2_list,
                        exclude_moieties=["alanine_sidechain"],
                        output_dir=output_dir,
                        output_prefix=f"{uploaded_file.name}_",
                    )

                    # Display results
                    st.success("Analysis complete!")

                    # Summary metrics
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Cofactor Atoms", len(cofactor_sphere))
                    with col2:
                        st.metric("PCS Atoms", len(pcs_atoms))
                    with col3:
                        st.metric("SCS Atoms", len(scs_atoms))

                    # Tabs for different views
                    tab1, tab2, tab3, tab4, tab5 = st.tabs(["3D Network", "Summary", "PCS Details", "SCS Details", "Downloads"])

                    with tab1:
                        st.subheader("Interactive 3D Coordination Network")
                        st.markdown("""
                        **Controls:**
                        - **Color by Coordination Sphere**: Cofactor (black), PCS (blue), SCS (fuchsia)
                        - **Color by Element**: C (black), N (blue), O (red), S (yellow), metals (various)
                        - **Backbone toggle**: Show/hide full residue backbone atoms
                        - **Background toggle**: Show/hide axes and grid
                        - **Hover** over atoms to see details including moiety classification
                        """)

                        # Generate the interactive 3D plot
                        fig = plot_interactive_modes_with_network(
                            structure=structure,
                            cofactor_atoms=cofactor_sphere,
                            pcs_atoms=pcs_atoms,
                            scs_atoms=scs_atoms,
                            bond_lookup_table=bond_lookup,
                            pdb_name=uploaded_file.name,
                            cofactor_resname=", ".join(cofactor_list),
                            atom_type_colors=atom_type_colors,
                            links_csv_path=None,  # No CSV file in webapp context
                            return_fig=True,
                        )
                        st.plotly_chart(fig, use_container_width=True)

                    with tab2:
                        st.subheader("Coordination Network Summary")

                        # Cofactor info
                        st.markdown("**Cofactor Atoms:**")
                        if cofactor_sphere:
                            cof_df = pd.DataFrame([
                                {
                                    "Residue": a["residue"],
                                    "Number": a["residue_number"],
                                    "Chain": a["chain"],
                                    "Atom": a["name"],
                                    "Element": a.get("element", ""),
                                }
                                for a in cofactor_sphere
                            ])
                            st.dataframe(cof_df, use_container_width=True, hide_index=True)
                        else:
                            st.warning("No cofactor atoms found with the specified residue names.")

                    with tab3:
                        st.subheader("Primary Coordination Sphere (PCS)")
                        if pcs_atoms:
                            pcs_df = pd.DataFrame([
                                {
                                    "Residue": a["residue"],
                                    "Number": a["residue_number"],
                                    "Chain": a["chain"],
                                    "Atom": a["name"],
                                    "Element": a.get("element", ""),
                                }
                                for a in pcs_atoms
                            ])

                            # Group by residue
                            st.markdown(f"**{len(pcs_atoms)} atoms across {pcs_df[['Residue', 'Number', 'Chain']].drop_duplicates().shape[0]} residues**")
                            st.dataframe(pcs_df, use_container_width=True, hide_index=True)

                            # Residue composition
                            st.markdown("**Residue Composition:**")
                            res_counts = pcs_df.groupby("Residue").size().sort_values(ascending=False)
                            st.bar_chart(res_counts)
                        else:
                            st.info("No PCS atoms found within the distance cutoff.")

                    with tab4:
                        st.subheader("Secondary Coordination Sphere (SCS)")
                        if scs_atoms:
                            scs_df = pd.DataFrame([
                                {
                                    "Residue": a["residue"],
                                    "Number": a["residue_number"],
                                    "Chain": a["chain"],
                                    "Atom": a["name"],
                                    "Element": a.get("element", ""),
                                }
                                for a in scs_atoms
                            ])

                            st.markdown(f"**{len(scs_atoms)} atoms across {scs_df[['Residue', 'Number', 'Chain']].drop_duplicates().shape[0]} residues**")
                            st.dataframe(scs_df, use_container_width=True, hide_index=True)

                            # Residue composition
                            st.markdown("**Residue Composition:**")
                            res_counts = scs_df.groupby("Residue").size().sort_values(ascending=False)
                            st.bar_chart(res_counts)
                        else:
                            st.info("No SCS atoms found.")

                    with tab5:
                        st.subheader("Download Results")

                        # Create downloadable CSVs
                        if cofactor_sphere:
                            cof_csv = pd.DataFrame([
                                {
                                    "residue": a["residue"],
                                    "residue_number": a["residue_number"],
                                    "chain": a["chain"],
                                    "atom": a["name"],
                                    "element": a.get("element", ""),
                                    "x": a["coordinates"][0],
                                    "y": a["coordinates"][1],
                                    "z": a["coordinates"][2],
                                }
                                for a in cofactor_sphere
                            ]).to_csv(index=False)
                            st.download_button(
                                "Download Cofactor Atoms (CSV)",
                                cof_csv,
                                file_name=f"{uploaded_file.name}_cofactor.csv",
                                mime="text/csv",
                            )

                        if pcs_atoms:
                            pcs_csv = pd.DataFrame([
                                {
                                    "residue": a["residue"],
                                    "residue_number": a["residue_number"],
                                    "chain": a["chain"],
                                    "atom": a["name"],
                                    "element": a.get("element", ""),
                                    "x": a["coordinates"][0],
                                    "y": a["coordinates"][1],
                                    "z": a["coordinates"][2],
                                }
                                for a in pcs_atoms
                            ]).to_csv(index=False)
                            st.download_button(
                                "Download PCS Atoms (CSV)",
                                pcs_csv,
                                file_name=f"{uploaded_file.name}_pcs.csv",
                                mime="text/csv",
                            )

                        if scs_atoms:
                            scs_csv = pd.DataFrame([
                                {
                                    "residue": a["residue"],
                                    "residue_number": a["residue_number"],
                                    "chain": a["chain"],
                                    "atom": a["name"],
                                    "element": a.get("element", ""),
                                    "x": a["coordinates"][0],
                                    "y": a["coordinates"][1],
                                    "z": a["coordinates"][2],
                                }
                                for a in scs_atoms
                            ]).to_csv(index=False)
                            st.download_button(
                                "Download SCS Atoms (CSV)",
                                scs_csv,
                                file_name=f"{uploaded_file.name}_scs.csv",
                                mime="text/csv",
                            )

                # Cleanup temp file
                os.unlink(tmp_path)

            except Exception as e:
                st.error(f"Error during analysis: {str(e)}")
                import traceback
                st.code(traceback.format_exc())

elif uploaded_file is None:
    # Show example/instructions when no file uploaded
    st.info("Upload a structure file and configure parameters in the sidebar to begin.")

    with st.expander("How to Use"):
        st.markdown("""
        ### Quick Start
        1. **Upload** a PDB or mmCIF structure file
        2. **Enter** the cofactor residue name(s) (e.g., `CU`, `HEM`, `ICS,CLF`)
        3. **Adjust** the distance cutoff if needed (default 3.6 Å)
        4. **Click** "Run Analysis"

        ### What This Tool Does
        This tool identifies the **coordination network** around protein cofactors:

        - **Primary Coordination Sphere (PCS)**: Atoms directly coordinating the cofactor
        - **Secondary Coordination Sphere (SCS)**: Atoms coordinating the PCS residues

        ### Example Cofactors
        | Cofactor | Residue Name | Example PDB |
        |----------|--------------|-------------|
        | Copper | CU | 1AG6 (Plastocyanin) |
        | Heme | HEM, HEA, HM1 | 1HHO (Hemoglobin) |
        | Iron-Sulfur | ICS, CLF | 3U7Q (Nitrogenase) |
        | OEC | OEX | 4UB6 (Photosystem II) |
        """)

    with st.expander("About SSCNA"):
        st.markdown("""
        **SSCNA** (Single Structure Cofactor Network Analysis) is part of the
        [Protein-Cofactor Coordination Network Analysis Tools](https://github.com/aaronledray/Protein-Cofactor_Coordination_Network_Analysis_Tools).

        It provides a robust, quantitative definition of coordination networks suitable for:
        - Protein structure analysis
        - Protein design workflows
        - Comparative analysis of metalloprotein active sites
        - AlphaFold3 model screening
        """)

# Footer
st.markdown("---")
st.markdown(
    "Built with [Streamlit](https://streamlit.io) | "
    "[GitHub](https://github.com/aaronledray/Protein-Cofactor_Coordination_Network_Analysis_Tools)"
)
