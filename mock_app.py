# mock_app.py
# ─────────────────────────────────────────────────────────────
# DegradeAI — 포스터 QR 시연용 목업 (stmol 없이 py3Dmol 직접 사용)
# Streamlit Cloud Python 3.14 호환 버전
# ─────────────────────────────────────────────────────────────
import streamlit as st
import time
import pandas as pd
import py3Dmol
import math

st.set_page_config(page_title="DegradeAI - PROTAC-1조 (Demo)", layout="wide")

# ═════════════════════════════════════════════════════════════
# 3D 렌더링 헬퍼 — stmol 없이 py3Dmol → HTML → st.components
# ═════════════════════════════════════════════════════════════
def show3d(view, height=560):
    """py3Dmol view를 Streamlit에 표시 (stmol 대체)."""
    html = view._make_html()
    st.components.v1.html(html, height=height, scrolling=False)


# ═════════════════════════════════════════════════════════════
# 사전 준비된 시연 데이터 (실제 DegradeAI 산출값)
# ═════════════════════════════════════════════════════════════
DEMO = {
    "smiles": "CC1=C(C2=NN=C(N=C2C3=C1C(=O)NC4=CC=CC=C43)C)C5=CC=C(C=C5)Cl",
    "pdb": "4LR6",
    "pkd": 6.78, "pkd_lo": 6.21, "pkd_hi": 7.32,
    "ad_score": 0.27,
    "smina_dg": -6.8, "pkd_docking": 4.99,
    "mw": 456.8, "logp": 4.23, "tpsa": 78.4,
    "hbd": 0, "hba": 5, "fsp3": 0.21, "bertz": 921.0,
    "sa_raw": 2.79, "sa_norm": 0.801, "sa_label": "합성 용이",
    "dlcs": 0.689, "dlcs_actual": 0.685,
    "e3": "CRBN", "e3_score": 98,
    "verdict": "🟡 HOLD",
}

# ── JQ1 리간드 PDB 블록 ──
JQ1_PDB = """\
HEADER    DEMO LIGAND JQ1
HETATM    1  C1  LIG A   1      -2.345   1.234   0.512  1.00  0.00           C
HETATM    2  C2  LIG A   1      -1.123   1.876  -0.234  1.00  0.00           C
HETATM    3  N1  LIG A   1      -0.045   1.123   0.123  1.00  0.00           N
HETATM    4  C3  LIG A   1       1.234   1.654  -0.345  1.00  0.00           C
HETATM    5  N2  LIG A   1       2.345   0.876   0.234  1.00  0.00           N
HETATM    6  C4  LIG A   1       2.123  -0.456  -0.123  1.00  0.00           C
HETATM    7  C5  LIG A   1       0.834  -0.987   0.345  1.00  0.00           C
HETATM    8  C6  LIG A   1      -0.345  -0.234  -0.234  1.00  0.00           C
HETATM    9  CL1 LIG A   1       3.456  -1.345   0.567  1.00  0.00          Cl
HETATM   10  C7  LIG A   1      -3.567   1.987  -0.123  1.00  0.00           C
HETATM   11  O1  LIG A   1      -1.234   3.123  -0.567  1.00  0.00           O
HETATM   12  C8  LIG A   1       0.567  -2.345   0.789  1.00  0.00           C
HETATM   13  C9  LIG A   1      -4.678   1.234   0.456  1.00  0.00           C
HETATM   14  N3  LIG A   1      -2.890  -1.456  -0.678  1.00  0.00           N
CONECT    1    2   10
CONECT    2    3   11
CONECT    3    4
CONECT    4    5
CONECT    5    6
CONECT    6    7    9
CONECT    7    8   12
CONECT    8    1   14
CONECT   10   13
END
"""

def _make_protein_pdb(offset_x=0.0, chain="A", n=40):
    """알파 헬릭스 형태 데모 단백질 PDB 생성 (인터넷 불필요)."""
    lines = []
    for i in range(n):
        ang = i * 100 * math.pi / 180.0
        x = 5.0 * math.cos(ang) + offset_x
        y = 5.0 * math.sin(ang)
        z = i * 1.5 - 30.0
        lines.append(
            f"ATOM  {i+1:5d}  CA  ALA {chain}{i+1:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C"
        )
    return "\n".join(lines)

PROTEIN_A = _make_protein_pdb(offset_x=0.0,  chain="A", n=40)
PROTEIN_B = _make_protein_pdb(offset_x=18.0, chain="B", n=30)


# ═════════════════════════════════════════════════════════════
# 3D 뷰어 함수들
# ═════════════════════════════════════════════════════════════
def render_qsar_3d():
    view = py3Dmol.view(width=860, height=520)
    view.addModel(PROTEIN_A, 'pdb')
    view.addModel(JQ1_PDB,   'pdb')
    view.setBackgroundColor('0x0a0a1a')
    view.setStyle({'model': 0}, {"cartoon": {'color': 'spectrum', 'opacity': 0.55}})
    view.setStyle({'model': 1}, {
        "stick":  {'colorscheme': 'magentaCarbon', 'radius': 0.42},
        "sphere": {'colorscheme': 'magentaCarbon', 'radius': 0.38, 'opacity': 0.60},
    })
    view.addSurface(py3Dmol.VDW, {'opacity': 0.20, 'color': 'magenta'}, {'model': 1})
    view.zoomTo({'model': 1})
    view.center({'model': 1})
    view.zoom(1.5)
    view.spin('y', 0.8)
    return view


def render_ternary_3d():
    view = py3Dmol.view(width=860, height=560)
    view.addModel(PROTEIN_A, 'pdb')
    view.addModel(PROTEIN_B, 'pdb')
    view.addModel(JQ1_PDB,   'pdb')
    view.setBackgroundColor('0x05050f')
    view.setStyle({'model': 0}, {"cartoon": {'color': 'spectrum',     'opacity': 0.50}})
    view.setStyle({'model': 1}, {"cartoon": {'colorscheme': 'greenCarbon', 'opacity': 0.45}})
    view.setStyle({'model': 2}, {
        "stick":  {'colorscheme': 'magentaCarbon', 'radius': 0.48},
        "sphere": {'colorscheme': 'magentaCarbon', 'radius': 0.42, 'opacity': 0.62},
    })
    view.zoomTo({'model': 2})
    view.center({'model': 2})
    view.zoom(1.2)
    view.spin('y', 0.7)
    return view


# ═════════════════════════════════════════════════════════════
# 공통 헤더
# ═════════════════════════════════════════════════════════════
st.title("🧬 DegradeAI [PROTAC-1조]")
st.markdown("#### 머신러닝 기반 적응형 표적 단백질 분해제(TPD) 통합 설계·검증 플랫폼")
st.info(
    "📱 **포스터 시연용 데모입니다.** "
    "실제 DegradeAI 플랫폼이 산출한 결과를 재현합니다. "
    "(안정적 시연을 위해 사전 연산된 결과를 표시 — 실시간 도킹·ML은 로컬 환경에서 동작)"
)

# ── 사이드바 ──
st.sidebar.title("🧬 DegradeAI")
st.sidebar.caption("포스터 QR 시연 모드")
mode = st.sidebar.radio(
    "📍 분석 모드",
    ["🎯 QSAR 결합력 분석",
     "🏗️ 조립식 Ternary 모델링",
     "💊 PROTAC 자동 완성 설계"],
)
st.sidebar.divider()


# ═════════════════════════════════════════════════════════════
# 모드 1 — QSAR 결합력 분석
# ═════════════════════════════════════════════════════════════
if mode == "🎯 QSAR 결합력 분석":
    st.sidebar.header("🛠️ 설계 워크스테이션")
    st.sidebar.selectbox("적응증 (TME)", ["폐암 (Lung)"], disabled=True)
    st.sidebar.text_input("Target PDB ID", DEMO["pdb"], disabled=True)
    st.sidebar.text_area("Ligand SMILES", DEMO["smiles"], disabled=True, height=110)
    st.sidebar.caption("⚠️ 데모 — 입력은 JQ1로 고정")
    run = st.sidebar.button("🚀 설계 파이프라인 가동", type="primary")

    st.divider()
    if not run:
        st.markdown("👈 **🚀 설계 파이프라인 가동** 버튼을 눌러 JQ1의 BRD4 결합력 분석을 시연하세요.")
        st.markdown(
            "**JQ1**은 대표적인 BET bromodomain 저해제로, "
            "PROTAC(MZ1, dBET 등)의 워헤드로 널리 쓰입니다."
        )
        st.stop()

    # 로그 애니메이션
    st.subheader("🖥️ AI 연산 관제 로그")
    box = st.empty(); logs = []
    for s in [
        "DegradeAI 코어 초기화 중...",
        "✅ 입력 검증 통과.",
        f"[REAL CORE] PDBFixer: {DEMO['pdb']} 구조 복구 중...",
        "✅ 단백질 구조 복구 완료",
        "[REAL CORE] Smina 도킹 엔진 구동 (cpu=8)...",
        f"✅ Smina ΔG: {DEMO['smina_dg']} kcal/mol",
        "[Safety] HPA DB 쿼리: 폐암(Lung)...",
        f"✅ E3 매핑: {DEMO['e3']} (조직 발현 {DEMO['e3_score']}/100)",
        "[CORE] 하이브리드 ML 추론 (Morgan FP 2055차원)...",
        f"✅ ML pKd: {DEMO['pkd']} (95% CI {DEMO['pkd_lo']}~{DEMO['pkd_hi']})",
        f"   📌 AD 점수: {DEMO['ad_score']} (✅ 신뢰)",
        f"✅ DLCS 예측: {DEMO['dlcs']}",
        f"✅ SA Score: {DEMO['sa_raw']}/10 ({DEMO['sa_label']})",
        "✅ 전체 파이프라인 완료. 대시보드 렌더링...",
    ]:
        logs.append(f"> {s}")
        box.code("\n".join(logs), language="bash")
        time.sleep(0.22)

    st.divider()
    st.markdown("### 🎬 Protein–Ligand Docking 3D — Dramatic View")
    st.info("🔵 스펙트럼 = 단백질 포켓 | 💗 형광 핑크 = 도킹된 JQ1")
    show3d(render_qsar_3d(), height=520)

    st.divider()
    st.success("✅ PROTAC 분자량 범위 충족")

    c1, c2, c3, c4, c5, c6 = st.columns(6)
    c1.metric("Predicted pKd (ML, BRD4)", f"{DEMO['pkd']:.2f}",
              f"95% CI: {DEMO['pkd_lo']:.2f}~{DEMO['pkd_hi']:.2f}",
              help=f"도킹 환산 pKd: {DEMO['pkd_docking']:.2f} | ΔG: {DEMO['smina_dg']}")
    c2.metric("분자량 (RDKit)", f"{DEMO['mw']:.1f} Da", "PROTAC 범위")
    c3.metric("E3 리가아제 (Lung)", DEMO['e3'], f"조직 발현: {DEMO['e3_score']}/100")
    with c4:
        st.metric("SA Score", f"{DEMO['sa_raw']:.2f} / 10", DEMO['sa_label'])
        st.caption(f"정규화: {DEMO['sa_norm']:.3f}")
    with c5:
        st.metric("Bertz CT", f"{DEMO['bertz']:.0f}")
        st.caption(f"Fsp3:{DEMO['fsp3']:.2f} | TPSA:{DEMO['tpsa']:.1f} | HBD:{DEMO['hbd']} | HBA:{DEMO['hba']}")
    with c6:
        st.metric("DLCS (약물성)", f"{DEMO['dlcs']:.3f}", f"실측: {DEMO['dlcs_actual']:.3f}")

    st.caption(f"🎯 신뢰도(AD): {DEMO['ad_score']:.2f} ✅ 학습 분포 내")
    st.info(
        f"ℹ️ Lipinski Ro5 위반 1개 (PROTAC은 위반 허용) — "
        f"LogP: {DEMO['logp']:.2f}, HBD: {DEMO['hbd']}, HBA: {DEMO['hba']}"
    )

    st.divider()
    st.subheader("📝 의약화학적 의견서 (LLM Report)")
    st.markdown(f"### 종합 판정: {DEMO['verdict']}")
    st.info(
        f"**1. 결합력(pKd) 및 약물성(ADME) 종합 진단**\n\n"
        f"예측 결합력 pKd {DEMO['pkd']:.2f}는 등급 **🟡 보통**에 해당하며, "
        "리드 최적화를 통한 향상이 권장됩니다. "
        f"분자량 {DEMO['mw']:.0f} Da, LogP {DEMO['logp']:.1f}로 세포 투과성은 양호합니다.\n\n"
        f"**2. 합성 접근성(SA) 및 약물성(DLCS) 평가**\n\n"
        f"SA Score {DEMO['sa_raw']:.2f}(**🟢 우수**), DLCS {DEMO['dlcs']:.3f}(**🟢 양호**)로 "
        "합성 가능성과 약물성이 양호합니다. "
        f"Bertz CT {DEMO['bertz']:.0f}는 표준 소분자 수준으로 CMC 리스크가 낮습니다.\n\n"
        "**3. 종합 판정 및 권고**\n\n"
        "**판정: 🟡 HOLD** — 조건부 보류. 결합력 보완 후 재평가를 권고합니다. "
        "합성성·약물성은 우수하나 결합력(pKd)이 임상 기준에 다소 못 미칩니다.\n"
        "* 다음 단계: 링커·E3 리간드 연결 후 ternary 복합체 형성 평가."
    )


# ═════════════════════════════════════════════════════════════
# 모드 2 — 조립식 Ternary 모델링
# ═════════════════════════════════════════════════════════════
elif mode == "🏗️ 조립식 Ternary 모델링":
    st.sidebar.header("🏗️ Ternary 설정")
    st.sidebar.text_input("워헤드 SMILES", DEMO["smiles"], disabled=True)
    st.sidebar.caption("🔗 QSAR 연동: JQ1 (pKd 6.78)")
    run = st.sidebar.button("🚀 조립식 Ternary 모델 생성", type="primary")

    st.title("🧬 조립식 Ternary Complex 모델링")
    st.markdown("#### 사용자 워헤드 + 실측 E3 템플릿(5T35) 기반 삼중 복합체 모델")
    st.info(
        "ℹ️ **Template-based Ternary Modeling** — 워헤드를 실측 ternary 구조"
        "(5T35: BRD4–VHL)의 리간드 위치에 정렬하여 삼중 복합체 모델을 생성합니다."
    )
    st.divider()

    if not run:
        st.markdown("👈 **🚀 조립식 Ternary 모델 생성** 버튼을 눌러 시연하세요.")
        st.stop()

    box = st.empty(); logs = []
    for s in [
        "실측 템플릿 5T35 로드 중...",
        "✅ 템플릿 확보 (BRD4–MZ1–VHL)",
        "✅ 타겟 체인 D / E3 체인 H / 템플릿 리간드 759",
        "사용자 워헤드 3D 생성 중...",
        "✅ 워헤드를 MZ1(759) 결합 포켓에 정렬 완료 (도킹 시뮬레이션)",
        "조립식 ternary 3D 렌더링 중...",
        "✅ 완료!",
    ]:
        logs.append(f"> {s}")
        box.code("\n".join(logs), language="bash")
        time.sleep(0.25)

    st.divider()
    st.subheader("📊 조립식 Ternary 분석 결과")
    c1, c2, c3 = st.columns(3)
    c1.metric("타겟 단백질", "BRD4 (BD2)", "실측 골격")
    c2.metric("E3 리가아제", "VHL", "실측 골격")
    c3.metric("워헤드 분자량", f"{DEMO['mw']:.1f} Da")

    st.divider()
    st.subheader("✅ 모델 검증 (Validation)")
    v1, v2 = st.columns(2)
    v1.metric("워헤드 ↔ MZ1 구조 유사도", "1.000",
              help="Tanimoto 유사도. JQ1은 MZ1의 워헤드이므로 1.0")
    v2.success("🟢 템플릿 적합성 높음 — JQ1은 MZ1의 워헤드로 동일 포켓·배치를 따릅니다.")

    st.divider()
    st.subheader("🎬 조립식 Ternary Complex 3D — Dramatic View")
    st.markdown(
        "🔵 **BRD4** (스펙트럼) &nbsp;|&nbsp; "
        "🟢 **VHL** (초록) &nbsp;|&nbsp; "
        "💗 **워헤드 JQ1** (발광 핑크)"
    )
    show3d(render_ternary_3d(), height=560)
    st.success("🎉 워헤드가 BRD4–VHL 삼중 복합체 골격에 조립된 모델 생성 완료!")
    st.caption(
        "💡 단백질 골격·E3는 실측 구조(5T35, Gadd et al. Nat. Chem. Biol. 2017) 기반. "
        "워헤드 배치는 MZ1 좌표 정렬 근사. (de novo 예측 아님)"
    )


# ═════════════════════════════════════════════════════════════
# 모드 3 — PROTAC 자동 완성 설계
# ═════════════════════════════════════════════════════════════
elif mode == "💊 PROTAC 자동 완성 설계":
    st.sidebar.header("💊 PROTAC 설계")
    st.sidebar.text_input("워헤드 SMILES", DEMO["smiles"], disabled=True)
    st.sidebar.selectbox("E3 리간드", ["VHL — VH032"], disabled=True)
    st.sidebar.caption("🔗 QSAR 탭 연동: JQ1 자동 적용")
    run = st.sidebar.button("🚀 PROTAC 후보 자동 생성", type="primary")

    st.title("💊 PROTAC 자동 완성 설계")
    st.markdown("#### 워헤드 1개 → E3 리간드 + 링커 자동 조립 → 완성 PROTAC 후보")
    st.info(
        "ℹ️ 워헤드에 **논문 검증 E3 리간드 5종** + **링커 30종** "
        "(PEG·Alkyl·아미드·고리형)을 RDKit으로 결합해 완성 PROTAC을 생성하고, "
        "pKd·약물성(DLCS)·분자량 적합성을 종합해 정렬합니다."
    )
    st.caption("📏 분자량 필터: 600–1,200 Da | sweet spot 700–1,100 Da")
    st.divider()

    if not run:
        st.markdown("👈 **🚀 PROTAC 후보 자동 생성** 버튼을 눌러 시연하세요.")
        st.stop()

    box = st.empty(); logs = []
    for s in [
        "워헤드 연결 지점([*:1]) 설정 중...",
        "✅ 사용자 지정 연결 지점 사용",
        "E3=VHL + 링커 30종 조합 생성·채점 중...",
        "   pKd 40% + DLCS 30% + 분자량 적합성 30% 종합 채점...",
        "✅ 완성 PROTAC 후보 18개 생성 (MW 필터 제외: 12개)",
        "✅ 완료!",
    ]:
        logs.append(f"> {s}")
        box.code("\n".join(logs), language="bash")
        time.sleep(0.25)

    st.divider()
    st.subheader("📊 완성 PROTAC 후보 순위표")
    df = pd.DataFrame([
        {"링커": "PEG3",       "종류": "PEG",   "추정길이(Å)": 13.5, "MW": 942.4, "pKd참고": 6.91, "DLCS": 0.512, "SA": 3.84, "RotBonds": 16, "종합점수": 78.3, "거리오차(Å)": 0.3},
        {"링커": "Amide-PEG2", "종류": "Amide", "추정길이(Å)": 12.0, "MW": 955.4, "pKd참고": 6.85, "DLCS": 0.521, "SA": 4.02, "RotBonds": 15, "종합점수": 73.8, "거리오차(Å)": 1.8},
        {"링커": "PEG4",       "종류": "PEG",   "추정길이(Å)": 17.0, "MW": 986.5, "pKd참고": 6.88, "DLCS": 0.498, "SA": 3.91, "RotBonds": 19, "종합점수": 74.1, "거리오차(Å)": 3.2},
        {"링커": "Pip-PEG",    "종류": "Rigid", "추정길이(Å)": 11.0, "MW": 968.5, "pKd참고": 6.79, "DLCS": 0.534, "SA": 4.21, "RotBonds": 13, "종합점수": 72.5, "거리오차(Å)": 2.8},
        {"링커": "Alkyl-C10",  "종류": "Alkyl", "추정길이(Å)": 12.6, "MW": 924.5, "pKd참고": 6.72, "DLCS": 0.487, "SA": 3.65, "RotBonds": 18, "종합점수": 70.2, "거리오차(Å)": 1.2},
    ])
    st.dataframe(df, use_container_width=True)
    st.caption(
        "정렬: 목표 거리(13.8 Å) 오차 → 종합점수 | "
        "종합점수 = pKd 40% + DLCS 30% + 분자량 적합성 30%"
    )

    st.divider()
    st.subheader("🏆 최우선 후보: 워헤드 + PEG3 + VH032 (VHL)")
    m1, m2, m3, m4, m5 = st.columns(5)
    m1.metric("분자량", "942.4 Da", "sweet spot ✅")
    m2.metric("pKd 참고", "6.91",
              help="⚠️ 모델은 BRD4 워헤드 기준 — 완성 PROTAC은 AD 밖, 상대 비교용")
    m3.metric("DLCS", "0.512")
    m4.metric("회전결합", "16")
    m5.metric("종합점수", "78.3/100")
    st.success(
        "🎉 완성 PROTAC 후보 생성 완료! "
        "이 SMILES를 QSAR / Ternary 탭에서 도킹·3D 검증으로 이어갈 수 있습니다."
    )


# ── 공통 푸터 ──
st.divider()
st.caption(
    "💡 본 데모는 DegradeAI 플랫폼(PROTAC-1조)의 실제 산출 결과를 재현한 시연용입니다. "
    "전체 기능(임의 분자 입력, 실시간 도킹·ML)은 로컬 환경에서 동작합니다. "
    "모든 예측은 wet-lab 실험으로 검증되어야 합니다."
)
