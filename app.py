import streamlit as st
import pandas as pd
import numpy as np
from utils.constants import HPO_FILE
from services.hpo_service import load_hpo_index, hpo_match
from services.vcf_service import parse_vcf
from services.vep_service import annotate_vep
from services.dbnsfp_service import (
    load_dbnsfp_annotations,
    merge_dbnsfp_annotations,
)
from services.ranking_service import rank_variants, score_label

st.set_page_config(
    page_title="Platforma pentru prioritizarea variantelor genetice",
    page_icon="🧬",
    layout="wide",
)

# ── CSS ───────────────────────────────────────────────────────────────────────
st.markdown("""
<style>
.stApp { max-width: 1200px; margin: 0 auto; }
.hpo-tag {
    display: inline-block;
    background: #E1F5EE; color: #085041;
    border: 1px solid #9FE1CB; border-radius: 6px;
    padding: 3px 10px; margin: 3px; font-size: 13px;
}
</style>
""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# INPUT PANEL
# ══════════════════════════════════════════════════════════════════════════════

def render_input_panel(hpo_terms_df):
    if "selected_hpo" not in st.session_state:
        st.session_state.selected_hpo = []

    if "vcf_content" not in st.session_state:
        st.session_state.vcf_content = None
    if "dbnsfp_basic_df" not in st.session_state:
        st.session_state.dbnsfp_basic_df = None
    if "dbnsfp_af_df" not in st.session_state:
        st.session_state.dbnsfp_af_df = None

    def build_dbnsfp_batches(vcf_content, min_qual, batch_size=20, max_variants=100):
        try:
            variants = parse_vcf(vcf_content, min_qual=min_qual).head(max_variants)
        except Exception:
            return []

        queries = [
            f"{str(row['chrom']).replace('chr', '')}-{row['pos']}-{row['ref']}-{row['alt']}"
            for _, row in variants.iterrows()
        ]
        return [
            "\n".join(queries[i:i + batch_size])
            for i in range(0, len(queries), batch_size)
        ]

    with st.container():
        st.subheader("1. Date pacient")

        c1, c2, c3 = st.columns(3)
        with c1:
            patient_id = st.text_input("ID pacient", value="Patient-001")
        with c2:
            patient_age = st.text_input("Varsta", value="20")
        with c3:
            patient_gender = st.selectbox("Sex", options=["Masculin", "Feminin"], index=0)

        genetician = st.text_input("Genetician", value="Dr. Ionescu")

        st.divider()

        st.subheader("2. Termeni HPO")
        st.caption("Cauta simptomele pacientului folosind termeni HPO.")

        st.markdown("**Scenarii predefinite:**")
        pcol1, pcol2 = st.columns(2)

        with pcol1:
            if st.button("Epilepsie pediatrica", use_container_width=True):
                st.session_state.selected_hpo = [
                    ("HP:0001250", "Seizure"),
                    ("HP:0000750", "Delayed speech and language development"),
                    ("HP:0001252", "Hypotonia"),
                ]
            if st.button("Distrofie musculara", use_container_width=True):
                st.session_state.selected_hpo = [
                    ("HP:0003560", "Muscular dystrophy"),
                    ("HP:0001252", "Hypotonia"),
                    ("HP:0008981", "Calf muscle hypertrophy"),
                ]
            if st.button("Retinita pigmentara", use_container_width=True):
                st.session_state.selected_hpo = [
                    ("HP:0001123", "Visual field defect"),
                    ("HP:0011463", "Childhood onset"),
                    ("HP:0000662", "Nyctalopia (NIGHT BLINDNESS)"),
                ]

        with pcol2:
            if st.button("Cardiomiopatie", use_container_width=True):
                st.session_state.selected_hpo = [
                    ("HP:0001639", "Hypertrophic cardiomyopathy"),
                    ("HP:0001671", "Abnormal cardiac septum morphology"),
                    ("HP:0004749", "Atrial flutter"),
                ]
            if st.button("Boala metabolica", use_container_width=True):
                st.session_state.selected_hpo = [
                    ("HP:0000822", "Hypertension"),
                    ("HP:0004322", "Short stature"),
                    ("HP:0000252", "Microcephaly"),
                ]

        st.markdown("**Sau cauta manual:**")
        search_q = st.text_input(
            "Cauta termen HPO",
            placeholder="ex: seizure, hypotonia, HP:0001250...",
            label_visibility="collapsed",
        )

        if search_q:
            q = search_q.lower().strip()
            mask = (
                hpo_terms_df["hpo_id"].str.lower().str.contains(q, na=False) |
                hpo_terms_df["hpo_name"].str.lower().str.contains(q, na=False)
            )
            matches = hpo_terms_df[mask].head(8)

            if not matches.empty:
                for _, row in matches.iterrows():
                    already = any(h[0] == row["hpo_id"] for h in st.session_state.selected_hpo)
                    label = f"{row['hpo_id']} — {row['hpo_name']} ({int(row['gene_count'])} gene)"
                    if already:
                        st.markdown(f"✅ ~~{label}~~")
                    else:
                        if st.button(f"➕ {label}", key=f"add_{row['hpo_id']}"):
                            st.session_state.selected_hpo.append((row["hpo_id"], row["hpo_name"]))
                            st.rerun()
            else:
                st.caption("Niciun rezultat gasit.")

        if st.session_state.selected_hpo:
            st.markdown("**Termeni selectati:**")
            tags_html = " ".join(
                f'<span class="hpo-tag"><b>{code}</b> {name}</span>'
                for code, name in st.session_state.selected_hpo
            )
            st.markdown(tags_html, unsafe_allow_html=True)

            to_remove = st.selectbox(
                "Sterge un termen:",
                options=["—"] + [f"{c} — {n}" for c, n in st.session_state.selected_hpo],
                label_visibility="collapsed",
            )
            if to_remove != "—":
                code_to_rm = to_remove.split(" — ")[0]
                st.session_state.selected_hpo = [
                    h for h in st.session_state.selected_hpo if h[0] != code_to_rm
                ]
                st.rerun()
        else:
            st.info("Adauga cel putin 1 termen HPO pentru a continua.")

        st.divider()

        st.subheader("3. Fisier VCF")
        st.caption("Incarca fisierul VCF al pacientului.")

        vcf_file = st.file_uploader("Incarca fisier VCF", type=["vcf", "txt"])

        if vcf_file:
            try:
                st.session_state.vcf_content = vcf_file.read().decode("utf-8")
                st.success(
                    f"{vcf_file.name} INCARCAT cu succes! "
                    f"({len(st.session_state.vcf_content.splitlines())} linii) ✓"
                )
            except Exception as e:
                st.error(f"Nu s-a putut incarca fisierul VCF: {e}")

        vcf_content = st.session_state.vcf_content

        st.divider()

        with st.expander("Setari avansate"):
            min_qual = st.slider("QUAL minim", 0, 60, 20)
            top_n = st.slider("Nr variante in top", 5, 100, 10)

        st.divider()

        st.subheader("4. Adnotari dbNSFP")

        if vcf_content:
            batches = build_dbnsfp_batches(vcf_content, min_qual=min_qual)
            if batches:
                st.markdown("**Query dbNSFP generat direct din VCF:**")
                st.caption(
                    "Foloseste aceste batch-uri in dbNSFP Web Query, descarca CSV-ul, apoi incarca-l mai jos. "
                    "Dupa asta apesi Incepe analiza o singura data."
                )
                for idx, batch in enumerate(batches, start=1):
                    st.text_area(
                        f"Batch dbNSFP {idx}/{len(batches)}",
                        value=batch,
                        height=120,
                        key=f"dbnsfp_pre_batch_{idx}",
                    )
                if len(batches) == 5:
                    st.caption("Sunt afisate primele 100 variante filtrate; pentru tot VCF-ul e mai potrivita baza dbNSFP locala.")

        dbnsfp_basic_file = st.file_uploader(
            "Incarca dbNSFP basic annotation",
            type=["tsv", "txt", "csv"],
            key="dbnsfp_basic_uploader",
        )

        dbnsfp_af_file = st.file_uploader(
            "Incarca dbNSFP allele frequency",
            type=["tsv", "txt", "csv"],
            key="dbnsfp_af_uploader",
        )

        if dbnsfp_basic_file is None:
            st.session_state.dbnsfp_basic_df = None
        else:
            try:
                st.session_state.dbnsfp_basic_df = load_dbnsfp_annotations(dbnsfp_basic_file)
                st.success(
                    f"Basic annotation incarcat: {dbnsfp_basic_file.name} "
                    f"({len(st.session_state.dbnsfp_basic_df)} variante)."
                )
            except Exception as e:
                st.session_state.dbnsfp_basic_df = None
                st.error(f"Nu s-a putut citi fisierul basic annotation: {e}")

        if dbnsfp_af_file is None:
            st.session_state.dbnsfp_af_df = None
        else:
            try:
                st.session_state.dbnsfp_af_df = load_dbnsfp_annotations(dbnsfp_af_file)
                popfreq_cols = [
                    c for c in st.session_state.dbnsfp_af_df.columns
                    if str(c).startswith("dbnsfp_popfreq_")
                ]
                st.success(
                    f"Allele frequency incarcat: {dbnsfp_af_file.name} "
                    f"({len(st.session_state.dbnsfp_af_df)} variante, {len(popfreq_cols)} coloane AF)."
                )
            except Exception as e:
                st.session_state.dbnsfp_af_df = None
                st.error(f"Nu s-a putut citi fisierul allele frequency: {e}")

        if st.session_state.dbnsfp_basic_df is None and st.session_state.dbnsfp_af_df is None:
            st.warning(f": Incarca un TSV/CSV exportat din dbNSFP pentru variantele pacientului") 
        else:
            parts = []
            if st.session_state.dbnsfp_basic_df is not None:
                parts.append("basic annotation")
            if st.session_state.dbnsfp_af_df is not None:
                parts.append("allele frequency")
            st.caption(f"dbNSFP activ: {', '.join(parts)} ")

        st.divider()

    return {
        "patient_id": patient_id,
        "patient_age": patient_age,
        "patient_gender": patient_gender,
        "genetician": genetician,
        "selected_hpo": st.session_state.selected_hpo,
        "vcf_content": vcf_content,
        "dbnsfp_basic_df": st.session_state.dbnsfp_basic_df,
        "dbnsfp_af_df": st.session_state.dbnsfp_af_df,
        "min_qual": min_qual,
        "top_n": top_n,
    }


def render_population_frequency(row):
    popfreq_cols = sorted([c for c in row.index if str(c).startswith("dbnsfp_popfreq_")])
    data = []
    seen_labels = set()
    for col in popfreq_cols:
        val = pd.to_numeric(pd.Series([row.get(col)]), errors="coerce").iloc[0]
        if pd.notna(val):
            label = str(col).replace("dbnsfp_popfreq_", "").replace("_", " ")
            seen_labels.add(label)
            data.append((label, float(val)))

    if pd.notna(row.get("dbnsfp_gnomad_af")) and "gnomad af" not in seen_labels:
        seen_labels.add("gnomad af")
        data.append(("gnomad af", float(row["dbnsfp_gnomad_af"])))

    if pd.notna(row.get("population_max_af")) and "max af" not in seen_labels:
        data.append(("max af", float(row["population_max_af"])))

    if not data:
        st.caption(
            "Aceasta varianta nu are frecvente populationale potrivite dupa chrom/pos/ref/alt "
            "in fisierul incarcat."
        )
        return

    freq_df = pd.DataFrame(data, columns=["Population", "AlleleFrequency"]).sort_values("AlleleFrequency", ascending=False)
    st.markdown("**Frecvență alelică pe populații**")
    st.bar_chart(freq_df.set_index("Population")["AlleleFrequency"])
    st.dataframe(
        freq_df.style.format({"AlleleFrequency": "{:.6f}"}),
        use_container_width=True,
        hide_index=True,
    )


# ══════════════════════════════════════════════════════════════════════════════
# LOAD DATA
# ══════════════════════════════════════════════════════════════════════════════

try:
    gene_to_hpo, hpo_terms_df, hpo_weight_map = load_hpo_index(HPO_FILE)
except Exception as e:
    st.error(str(e))
    st.stop()


# ══════════════════════════════════════════════════════════════════════════════
# LAYOUT
# ══════════════════════════════════════════════════════════════════════════════

left, right = st.columns([1, 1.4], gap="large")

with left:
    ui_data = render_input_panel(hpo_terms_df)

with right:
    selected_hpo = ui_data["selected_hpo"]
    vcf_content = ui_data["vcf_content"]
    dbnsfp_basic_df = ui_data["dbnsfp_basic_df"]
    dbnsfp_af_df = ui_data["dbnsfp_af_df"]
    min_qual = ui_data["min_qual"]
    top_n = ui_data["top_n"]

    can_run = (
        len(selected_hpo) >= 1 and
        vcf_content is not None
    )

    run_clicked = st.button(
        "Incepe analiza",
        type="primary",
        use_container_width=True,
        disabled=not can_run,
    )

    if not run_clicked and "last_results" not in st.session_state:
        st.markdown("""
        <div style="text-align:center; padding: 60px 20px; color: #8C8880;">
            <div style="font-size:16px; margin-bottom:8px">
                Completeaza datele si apasa <b>Incepe analiza</b>
            </div>
        </div>
        """, unsafe_allow_html=True)

    if run_clicked:
        try:
            patient_hpo_codes = [code for code, _ in selected_hpo]

            progress_bar = st.progress(0)
            status_text = st.empty()

            def update_progress(val, msg):
                progress_bar.progress(min(max(val, 0.0), 1.0))
                status_text.markdown(f"⏳ {msg}")

            update_progress(0.05, "Pasul 1/4 — Parsare VCF...")
            variants_df = parse_vcf(vcf_content, min_qual=min_qual)

            if variants_df.empty:
                raise ValueError("VCF-ul nu contine variante valide dupa filtrare (QUAL/FILTER).")

            update_progress(0.20, f"Pasul 1/4 — {len(variants_df)} variante bune")

            update_progress(0.28, "Pasul 2/4 — Incepe procesul de adnotare...")
            annotated_df = annotate_vep(
                variants_df,
                progress_cb=lambda v, m: update_progress(0.28 + v * 0.25, f"Pasul 2/4 — {m}")
            )

            update_progress(0.58, f"Pasul 2/4 — {len(annotated_df)} variante adnotate")

            update_progress(0.68, f"Pasul 3/4 — Potrivire HPO ({len(patient_hpo_codes)} termeni)...")
            matched_df = hpo_match(annotated_df, patient_hpo_codes, gene_to_hpo, hpo_weight_map)
            n_matched = int(matched_df["hpo_match"].sum())
            update_progress(0.82, f"Pasul 3/4 - {n_matched}/{len(matched_df)} variante cu potrivire HPO")

            update_progress(0.86, "Pasul 4/5 — Integrare dbNSFP...")
            annotated_with_dbnsfp = matched_df
            if dbnsfp_basic_df is not None and not dbnsfp_basic_df.empty:
                annotated_with_dbnsfp = merge_dbnsfp_annotations(annotated_with_dbnsfp, dbnsfp_basic_df)
            if dbnsfp_af_df is not None and not dbnsfp_af_df.empty:
                annotated_with_dbnsfp = merge_dbnsfp_annotations(annotated_with_dbnsfp, dbnsfp_af_df)

            dbnsfp_cols = [c for c in annotated_with_dbnsfp.columns if c.startswith("dbnsfp_")]
            dbnsfp_hits = int(annotated_with_dbnsfp[dbnsfp_cols].notna().any(axis=1).sum()) if dbnsfp_cols else 0
            popfreq_cols = [c for c in annotated_with_dbnsfp.columns if c.startswith("dbnsfp_popfreq_")]
            popfreq_hits = int(annotated_with_dbnsfp[popfreq_cols].notna().any(axis=1).sum()) if popfreq_cols else 0

            update_progress(0.92, f"Pasul 5/5 — Ranking variante (dbNSFP: {dbnsfp_hits}, AF: {popfreq_hits})...")
            results_df = rank_variants(annotated_with_dbnsfp)

            update_progress(1.0, "Analiza finalizata!")

            progress_bar.empty()
            status_text.empty()

            st.session_state.last_results = results_df
            st.session_state.last_patient = {
                "id": ui_data["patient_id"],
                "age": ui_data["patient_age"],
                "gender": ui_data["patient_gender"],
                "genetician": ui_data["genetician"],
                "hpo": selected_hpo,
                "n_variants": len(variants_df),
                "dbnsfp_hits": dbnsfp_hits,
                "popfreq_hits": popfreq_hits,
            }

        except Exception as e:
            st.error(f"Eroare în analiza: {e}")

    if "last_results" in st.session_state:
        results_df = st.session_state.last_results
        pat = st.session_state.last_patient
        top_df = results_df.head(top_n)

        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Pacient", pat["id"])
        m2.metric("Variante analizate", pat["n_variants"])
        m3.metric("Termeni HPO", len(pat["hpo"]))
        m4.metric("HPO matches", int(results_df["hpo_match"].sum()))
        st.caption(
            f"dbNSFP matches: {pat.get('dbnsfp_hits', 0)} | "
            f"variante cu frecvente populationale: {pat.get('popfreq_hits', 0)}"
        )

        st.divider()
        st.markdown(f"#### Top {len(top_df)} variante prioritizate")

        for i, row in top_df.iterrows():
            score = float(row["pathogenicity_score"])
            rank = i + 1
            score_icon, _ = score_label(score)

            hpo_tags = ""
            if row.get("matched_hpo_terms"):
                for term in str(row["matched_hpo_terms"]).split("|"):
                    if term:
                        hpo_tags += f'<span class="hpo-tag">{term}</span> '

            gene_symbol = row.get("gene_symbol")
            consequence = row.get("consequence")
            display_gene = gene_symbol if pd.notna(gene_symbol) and gene_symbol else "necunoscut"
            display_cons = consequence.replace("_", " ") if pd.notna(consequence) and consequence else "necunoscut"

            with st.expander(
                f"#{rank} **{display_gene}** {score_icon} scor: **{score:.4f}** — {display_cons}",
                expanded=(rank <= 3),
            ):
                ca, cb = st.columns(2)

                with ca:
                    st.markdown("**Pozitie genomica**")
                    st.code(
                        f"{row.get('chrom', '?')}:{row.get('pos', '?')}  {row.get('ref', '?')} → {row.get('alt', '?')}",
                        language="text"
                    )

                    if pd.notna(row.get("hgvsp")) and row.get("hgvsp"):
                        st.markdown(f"**Schimbare proteica:** `{row['hgvsp']}`")

                    if pd.notna(row.get("impact")) and row.get("impact"):
                        st.markdown(f"**Impact:** `{row['impact']}`")

                    if pd.notna(row.get("am_score")):
                        st.markdown(
                            f"**AlphaMissense:** `{float(row['am_score']):.3f}` — "
                            f"*{row['am_class'] if pd.notna(row.get('am_class')) else '?'}*"
                        )

                    if pd.notna(row.get("cadd_from_vep")):
                        st.markdown(f"**CADD din VEP:** `{float(row['cadd_from_vep']):.2f}`")

                    if pd.notna(row.get("dbnsfp_revel_score")):
                        st.markdown(f"**REVEL dbNSFP:** `{float(row['dbnsfp_revel_score']):.3f}`")

                    if pd.notna(row.get("dbnsfp_cadd_phred")):
                        st.markdown(f"**CADD dbNSFP:** `{float(row['dbnsfp_cadd_phred']):.2f}`")

                    if pd.notna(row.get("dbnsfp_alphamissense_score")):
                        st.markdown(f"**AlphaMissense dbNSFP:** `{float(row['dbnsfp_alphamissense_score']):.3f}`")

                    if pd.notna(row.get("dbnsfp_metarnn_score")):
                        st.markdown(f"**MetaRNN dbNSFP:** `{float(row['dbnsfp_metarnn_score']):.3f}`")

                    if pd.notna(row.get("dbnsfp_gnomad_af")):
                        st.markdown(f"**gnomAD AF dbNSFP:** `{float(row['dbnsfp_gnomad_af']):.5f}`")

                    if pd.notna(row.get("dbnsfp_gerp_rs")):
                        st.markdown(f"**GERP++ dbNSFP:** `{float(row['dbnsfp_gerp_rs']):.2f}`")

                    if pd.notna(row.get("population_max_af")):
                        st.markdown(f"**Max AF populațional:** `{float(row['population_max_af']):.6f}`")

                with cb:
                    st.markdown("**Factori de ranking**")
                    st.markdown(f"HPO match count: **{int(row.get('hpo_match_count', 0))}**")
                    if pd.notna(row.get("phenotype_score")):
                        st.markdown(f"Phenotype score: **{float(row['phenotype_score']):.3f}**")
                    if pd.notna(row.get("phenotype_block")):
                        st.markdown(f"Phenotype block: **{float(row['phenotype_block']):.3f}**")
                    if pd.notna(row.get("molecular_block")):
                        st.markdown(f"Molecular block: **{float(row['molecular_block']):.3f}**")
                    if pd.notna(row.get("consequence_block")):
                        st.markdown(f"Consequence block: **{float(row['consequence_block']):.3f}**")
                    if pd.notna(row.get("population_block")):
                        st.markdown(f"Population block: **{float(row['population_block']):.3f}**")

                    cons_score = pd.to_numeric(pd.Series([row.get("consequence_score")]), errors="coerce").iloc[0]
                    if pd.notna(cons_score):
                        st.markdown(f"Consequence score: **{float(cons_score):.2f}**")

                    impact_score = pd.to_numeric(pd.Series([row.get('impact_score')]), errors="coerce").iloc[0]
                    if pd.notna(impact_score):
                        st.markdown(f"Impact score: **{float(impact_score):.2f}**")

                    if pd.notna(row.get("revel_norm")):
                        st.markdown(f"REVEL norm: **{float(row['revel_norm']):.2f}**")

                    if pd.notna(row.get("sift_norm")):
                        st.markdown(f"SIFT norm: **{float(row['sift_norm']):.2f}**")

                    if pd.notna(row.get("polyphen_norm")):
                        st.markdown(f"PolyPhen norm: **{float(row['polyphen_norm']):.2f}**")

                    if pd.notna(row.get("metarnn_norm")):
                        st.markdown(f"MetaRNN norm: **{float(row['metarnn_norm']):.2f}**")

                if hpo_tags:
                    st.markdown(
                        f"**Termeni HPO potriviti ({int(row.get('hpo_match_count', 0))}):** {hpo_tags}",
                        unsafe_allow_html=True
                    )
                else:
                    st.markdown("**Termeni HPO potriviti:** —")

                st.progress(min(max(score, 0.0), 1.0))
                st.caption(f"Scor final ranking: {score:.4f}")

                st.divider()
                render_population_frequency(row)

        st.divider()

        st.markdown("#### Export rezultate")
        export_cols = [
            "gene_symbol", "chrom", "pos", "ref", "alt",
            "consequence", "impact", "hgvsp",
            "am_score", "am_class", "cadd_from_vep",
            "dbnsfp_revel_score", "dbnsfp_cadd_phred",
            "dbnsfp_alphamissense_score", "dbnsfp_metarnn_score",
            "dbnsfp_sift_score",
            "dbnsfp_polyphen2_hdiv_score", "dbnsfp_gerp_rs",
            "dbnsfp_gnomad_af", "population_max_af", "dbnsfp_aapos", "dbnsfp_hgvsp",
            "phenotype_score", "phenotype_block", "molecular_block",
            "consequence_block", "population_block",
            "hpo_match_count", "matched_hpo_terms",
            "pathogenicity_score"
        ]
        export_cols.extend([c for c in top_df.columns if str(c).startswith("dbnsfp_popfreq_")])

        export_cols = [c for c in export_cols if c in top_df.columns]
        csv = top_df[export_cols].to_csv(index=False)

        st.download_button(
            "Descarca Top rezultate (CSV)",
            data=csv,
            file_name=f"results_{pat['id']}.csv",
            mime="text/csv",
            use_container_width=True,
        )
