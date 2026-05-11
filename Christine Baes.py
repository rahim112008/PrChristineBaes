import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from io import BytesIO
from PIL import Image
import base64
import time

# Configuration de la page
st.set_page_config(page_title="LiveGene Suite", page_icon="🧬", layout="wide")

# ----- Barre latérale -----
st.sidebar.title("LiveGene Suite")
st.sidebar.markdown("*Pre Christine Baes - Génomique du bétail*")
module = st.sidebar.radio(
    "Modules",
    ["🏠 PhenoCollect", "🧪 GenoPipeline", "📊 SelectSim", "🧬 ConsangWatch", "📸 Morphométrie"]
)

# ----- Utilitaires -----
def telecharger_csv(df, nom_fichier):
    csv = df.to_csv(index=False).encode()
    st.download_button("📥 Télécharger CSV", csv, nom_fichier, "text/csv")

def telecharger_rapport(contenu, nom_fichier):
    st.download_button("📄 Télécharger rapport", contenu, nom_fichier, "text/plain")

# ----- Sessions -----
if 'pheno_data' not in st.session_state:
    st.session_state.pheno_data = pd.DataFrame(columns=["Date", "Animal", "Phénotype", "Valeur brute", "Insémination", "Valeur corrigée", "Analyse"])

if 'mesures_morpho' not in st.session_state:
    st.session_state.mesures_morpho = pd.DataFrame(columns=["Date", "Animal", "Mesure #", "Type", "Distance (px)", "Distance (cm)"])

# ========== PHENOCOLLECT ==========
if module == "🏠 PhenoCollect":
    st.header("📸 PhenoCollect – Saisie et correction des phénotypes")
    with st.form("pheno_form"):
        col1, col2, col3 = st.columns(3)
        with col1:
            animal_id = st.text_input("ID Animal", placeholder="HOLCAN0001")
            pheno_type = st.selectbox("Type de phénotype",
                ["distance_anogenitale", "taux_conception", "score_boiterie", "etat_corporel", "position_uterus"])
        with col2:
            pheno_value = st.text_input("Valeur / Fichier", placeholder="12.5")
            if pheno_type in ["distance_anogenitale", "taux_conception"]:
                insem_type = st.selectbox("Type d'insémination", ["naturelle", "IA", "TAI"])
            else:
                insem_type = "—"
                st.selectbox("Type d'insémination", ["—"], disabled=True)

        submitted = st.form_submit_button("🔍 Enregistrer & Analyser")

        if submitted and animal_id and pheno_value:
            raw_val = pheno_value
            try:
                num_val = float(pheno_value)
            except:
                num_val = None

            corrige = raw_val
            analyse = ""
            correction_appliquee = False

            if num_val is not None and pheno_type in ["distance_anogenitale", "taux_conception"]:
                if insem_type == "TAI":
                    correction_appliquee = True
                    if pheno_type == "distance_anogenitale":
                        corrige = f"{(num_val * 1.15):.2f}"
                    elif pheno_type == "taux_conception":
                        corrige = f"{(num_val + 10):.1f}"

            if num_val is not None:
                if pheno_type == "distance_anogenitale":
                    analyse = "Fertilité élevée" if num_val > 10 else "Fertilité modérée"
                elif pheno_type == "taux_conception":
                    analyse = "Bonne fertilité" if num_val >= 40 else "Fertilité à améliorer"
                elif pheno_type == "score_boiterie":
                    analyse = "Bonne locomotion" if num_val <= 2 else "Boiterie suspectée"
                else:
                    analyse = "Indicateur enregistré"
            else:
                analyse = "Image analysée (simulation)"

            if correction_appliquee:
                analyse += f" | ⚠️ Correction TAI : {raw_val} → {corrige}"

            new_entry = pd.DataFrame([[pd.Timestamp.now().strftime("%Y-%m-%d %H:%M"), animal_id, pheno_type, raw_val,
                                       insem_type if pheno_type in ["distance_anogenitale", "taux_conception"] else "—",
                                       corrige, analyse]],
                                     columns=st.session_state.pheno_data.columns)
            st.session_state.pheno_data = pd.concat([st.session_state.pheno_data, new_entry], ignore_index=True)
            st.success(f"Phénotype ajouté pour {animal_id}. {analyse}")

    st.subheader("Historique des entrées")
    if not st.session_state.pheno_data.empty:
        st.dataframe(st.session_state.pheno_data)
        telecharger_csv(st.session_state.pheno_data, "phenotypes_corriges.csv")
    else:
        st.info("Aucune entrée.")

# ========== GENOPIPELINE ==========
elif module == "🧪 GenoPipeline":
    st.header("🧪 GenoPipeline – Pipeline génomique automatique")
    st.markdown("Simulez l'imputation et une GWAS.")
    geno_file = st.file_uploader("Fichier PLINK/VCF (facultatif)", type=["ped", "map", "vcf"])
    if st.button("🚀 Lancer le pipeline"):
        progress_bar = st.progress(0)
        status_text = st.empty()
        log_console = st.empty()
        logs = ""

        etapes = [
            ("1/5 - Contrôle qualité (call rate, MAF, HWE)...", 10, "QC : 150 000 variants, 2 000 individus."),
            ("2/5 - Phasing avec Eagle2.4...", 25, "Phasing terminé, 5 000 haplotypes."),
            ("3/5 - Imputation FImpute vers séquence...", 50, "Imputation : 25M variants, concordance 0.97."),
            ("4/5 - Post-imputation QC (info > 0.4)...", 75, "18M variants retenus."),
            ("5/5 - GWAS (ssGBLUP) & prédictions...", 95, "Top SNP identifié sur fertilité."),
            ("Pipeline terminé.", 100, "Rapport disponible.")
        ]

        for msg, pct, log in etapes:
            status_text.text(msg)
            progress_bar.progress(pct)
            logs += f"[{pd.Timestamp.now().strftime('%H:%M:%S')}] {log}\n"
            log_console.code(logs)
            time.sleep(0.8)

        top_snp = f"rs{np.random.randint(10**7, 9*10**7)}"
        st.success("✅ Pipeline terminé")
        st.write(f"**Top SNP détecté :** {top_snp} (associé à la fertilité, p = 1.2e-8)")

        fig, ax = plt.subplots(figsize=(10, 4))
        chr_labels = [f"chr{i}" for i in range(1, 31)]
        pvals = -np.log10(np.random.rand(30))
        colors = ['#e74c3c' if p > 2.5 else '#27ae60' for p in pvals]
        ax.bar(chr_labels, pvals, color=colors)
        ax.set_ylabel("-log10(p)")
        ax.set_title("Manhattan Plot simulé")
        plt.xticks(rotation=45)
        st.pyplot(fig)

        rapport = f"RAPPORT D'ANALYSE GENOMIQUE\nDate : {pd.Timestamp.now()}\nPipeline : Eagle + FImpute + GWAS\nTop SNP : {top_snp}\nConclusion : variant candidat pour la fertilité."
        telecharger_rapport(rapport, "rapport_genomique.txt")

# ========== SELECTSIM ==========
elif module == "📊 SelectSim":
    st.header("📊 SelectSim – Simulateur de sélection durable")
    col1, col2 = st.columns([2, 1])
    with col1:
        w_lait = st.slider("Production laitière (%)", 0, 100, 50)
        w_fert = st.slider("Fertilité (%)", 0, 100, 30)
        w_sante = st.slider("Santé & Bien-être (%)", 0, 100 - w_lait - w_fert, 20)
        if st.button("📈 Simuler"):
            gain_lait = w_lait * 8
            delta_fert = (w_fert * 0.05) - 1.0
            consang = 0.02 + (100 - w_lait) * 0.0003
            with col2:
                st.metric("Gain lait (kg/an)", f"+{gain_lait:.0f}")
                st.metric("Évolution fertilité", f"{delta_fert:.2f} unités")
                st.metric("Consanguinité (F)", f"{consang:.3f}")
            rapport_sim = f"SIMULATION DE SÉLECTION\nDate : {pd.Timestamp.now()}\nPondérations : Lait {w_lait}%, Fertilité {w_fert}%, Santé {w_sante}%\nRésultats à 20 ans :\n- Gain lait : +{gain_lait:.0f} kg\n- Fertilité : {delta_fert:.2f}\n- Consanguinité : {consang:.3f}"
            telecharger_rapport(rapport_sim, "simulation_selection.txt")

# ========== CONSANGWATCH ==========
elif module == "🧬 ConsangWatch":
    st.header("🧬 ConsangWatch – Monitoring de la consanguinité")
    st.subheader("Tendance F_ROH")
    annees = [2018, 2019, 2020, 2021, 2022, 2023]
    f_roh = [0.015, 0.018, 0.022, 0.025, 0.029, 0.033]
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=annees, y=f_roh, mode='lines+markers', name='F_ROH moyen'))
    fig.update_layout(yaxis_title="F_ROH", xaxis_title="Année", yaxis_range=[0, 0.04])
    st.plotly_chart(fig, use_container_width=True)
    telecharger_csv(pd.DataFrame({"Année": annees, "F_ROH": f_roh}), "consanguinite.csv")

    st.subheader("Planification d'accouplements")
    pere = st.text_input("ID Père", "TAUREAU001")
    mere = st.text_input("ID Mère", "VACHE123")
    if st.button("🔎 Vérifier consanguinité prévue"):
        f_pred = np.random.uniform(0, 0.08)
        if f_pred > 0.05:
            st.error(f"⚠️ Accouplement à risque ! F_ROH prévu = {f_pred:.3f}")
        else:
            st.success(f"✅ Accouplement acceptable. F_ROH prévu = {f_pred:.3f}")

# ========== MORPHOMÉTRIE (HTML/JS autonome – sans drawable-canvas) ==========
elif module == "📸 Morphométrie":
    st.header("📸 Morphométrie – Mesures morphométriques sur animal")
    st.markdown("""
    **Étalonnez** l'échelle, choisissez votre mode de mesure, puis cliquez sur l'image.
    - **Droite** : 2 points → distance rectiligne.
    - **Courbe** : plusieurs points → somme des segments (idéal pour suivre le dos, le thorax…).
    La distance s'affiche directement sur l'image. Ensuite, recopiez la valeur dans le champ ci‑dessous pour l'enregistrer.
    """)

    animal_morpho = st.text_input("🐄 ID Animal", value="", placeholder="Ex: HOLCAN1234")

    mode = st.radio("Mode d'acquisition", ["📷 Prendre une photo", "📁 Télécharger une image"], index=0)
    image = None
    if mode == "📁 Télécharger une image":
        uploaded_file = st.file_uploader("Choisissez une image", type=["jpg", "jpeg", "png"])
        if uploaded_file is not None:
            image = Image.open(uploaded_file)
    else:
        camera_file = st.camera_input("Prenez une photo")
        if camera_file is not None:
            image = Image.open(camera_file)

    if image is not None:
        # Étalonnage
        st.subheader("⚖️ Étalonnage")
        col_ref1, col_ref2 = st.columns(2)
        with col_ref1:
            ref_px = st.number_input("Distance de référence (pixels)", min_value=1.0, value=100.0, step=1.0, key="ref_px")
        with col_ref2:
            ref_cm = st.number_input("Distance réelle correspondante (cm)", min_value=0.1, value=10.0, step=0.1, key="ref_cm")

        # Choix du type de mesure
        mesure_type = st.radio("Type de mesure", ["📏 Droite (2 points)", "〰️ Courbe (plusieurs points)"], index=0)
        max_points = 2 if mesure_type.startswith("📏") else 20  # Nombre max de points pour le mode courbe

        # Conversion de l'image en base64
        buffered = BytesIO()
        image.save(buffered, format="PNG")
        img_b64 = base64.b64encode(buffered.getvalue()).decode()
        img_w, img_h = image.size

        # HTML / JS du canevas
        html_canvas = f"""
        <div style="text-align:center;">
            <canvas id="morphoCanvas" width="{img_w}" height="{img_h}" style="border:1px solid #ccc; max-width:100%; height:auto; cursor:crosshair;"></canvas>
            <br>
            <button onclick="resetPoints()" style="margin:10px;">🗑️ Effacer les points</button>
            <p id="distanceDisplay" style="font-weight:bold;"></p>
        </div>
        <script>
            const canvas = document.getElementById('morphoCanvas');
            const ctx = canvas.getContext('2d');
            const img = new Image();
            img.onload = function() {{
                ctx.drawImage(img, 0, 0);
            }};
            img.src = "data:image/png;base64,{img_b64}";

            let points = [];
            const MAX_POINTS = {max_points};
            const refPx = {ref_px};
            const refCm = {ref_cm};

            canvas.addEventListener('click', function(e) {{
                const rect = canvas.getBoundingClientRect();
                const scaleX = canvas.width / rect.width;
                const scaleY = canvas.height / rect.height;
                const x = (e.clientX - rect.left) * scaleX;
                const y = (e.clientY - rect.top) * scaleY;
                points.push({{x, y}});
                if (points.length > MAX_POINTS) points = points.slice(-MAX_POINTS);
                redraw();
                updateDistance();
            }});

            function redraw() {{
                ctx.clearRect(0, 0, canvas.width, canvas.height);
                ctx.drawImage(img, 0, 0);
                points.forEach((p, i) => {{
                    ctx.fillStyle = 'red';
                    ctx.beginPath();
                    ctx.arc(p.x, p.y, 4, 0, 2*Math.PI);
                    ctx.fill();
                    ctx.fillStyle = 'white';
                    ctx.font = '12px Arial';
                    ctx.fillText('P'+(i+1), p.x+6, p.y-6);
                }});
                if (points.length >= 2) {{
                    ctx.strokeStyle = 'lime';
                    ctx.lineWidth = 2;
                    for (let i=0; i<points.length-1; i++) {{
                        ctx.beginPath();
                        ctx.moveTo(points[i].x, points[i].y);
                        ctx.lineTo(points[i+1].x, points[i+1].y);
                        ctx.stroke();
                    }}
                }}
            }}

            function updateDistance() {{
                if (points.length >= 2) {{
                    let totalPx = 0;
                    for (let i=0; i<points.length-1; i++) {{
                        const dx = points[i+1].x - points[i].x;
                        const dy = points[i+1].y - points[i].y;
                        totalPx += Math.sqrt(dx*dx + dy*dy);
                    }}
                    const distCm = (totalPx * refCm / refPx).toFixed(2);
                    document.getElementById('distanceDisplay').innerHTML = 
                        `Distance mesurée : ${{totalPx.toFixed(1)}} pixels → ${{distCm}} cm`;
                }} else {{
                    document.getElementById('distanceDisplay').innerHTML = '';
                }}
            }}

            function resetPoints() {{
                points = [];
                redraw();
                updateDistance();
            }}
        </script>
        """

        st.components.v1.html(html_canvas, height=img_h + 120, scrolling=False)

        # Champ manuel pour récupérer la mesure (en pixels) et l'enregistrer
        st.markdown("---")
        st.subheader("💾 Enregistrer la mesure")
        distance_px_input = st.number_input("Distance mesurée en pixels (recopiez la valeur affichée ci‑dessus)", min_value=0.0, step=0.1, format="%.1f")
        distance_cm_calc = (distance_px_input * ref_cm) / ref_px if ref_px != 0 else 0.0
        st.caption(f"Conversion automatique : {distance_px_input:.1f} px → {distance_cm_calc:.2f} cm")

        if st.button("💾 Enregistrer cette mesure"):
            if animal_morpho.strip() == "":
                st.warning("⚠️ Veuillez entrer un ID animal.")
            elif distance_px_input <= 0:
                st.warning("⚠️ Veuillez entrer une distance supérieure à 0.")
            else:
                type_label = "droite" if max_points == 2 else "courbe"
                n_mesure = len(st.session_state.mesures_morpho[st.session_state.mesures_morpho["Animal"] == animal_morpho]) + 1
                new_mesure = pd.DataFrame([[
                    pd.Timestamp.now().strftime("%Y-%m-%d %H:%M"),
                    animal_morpho,
                    n_mesure,
                    type_label,
                    round(distance_px_input, 2),
                    round(distance_cm_calc, 2)
                ]], columns=st.session_state.mesures_morpho.columns)
                st.session_state.mesures_morpho = pd.concat([st.session_state.mesures_morpho, new_mesure], ignore_index=True)
                st.success(f"Mesure #{n_mesure} enregistrée pour {animal_morpho}.")
                st.rerun()

    # --- Historique des mesures ---
    st.subheader("📋 Mesures enregistrées")
    if not st.session_state.mesures_morpho.empty:
        animal_filtre = st.text_input("🔍 Filtrer par ID Animal", "")
        if animal_filtre:
            df_aff = st.session_state.mesures_morpho[st.session_state.mesures_morpho["Animal"].str.contains(animal_filtre)]
        else:
            df_aff = st.session_state.mesures_morpho
        st.dataframe(df_aff, use_container_width=True)
        telecharger_csv(df_aff, "mesures_morpho.csv")
        if st.button("🗑️ Supprimer la dernière mesure enregistrée"):
            if not st.session_state.mesures_morpho.empty:
                st.session_state.mesures_morpho = st.session_state.mesures_morpho.iloc[:-1]
                st.rerun()
    else:
        st.info("Aucune mesure enregistrée.")
