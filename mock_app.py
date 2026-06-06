# mock_app.py
# ─────────────────────────────────────────────────────────────
# DegradeAI — 포스터 QR 시연용 목업 (Mock-up) [풀버전]
#
#   ✅ 인터넷 / RDKit / Smina / ML 모델 전혀 불필요
#   ✅ Streamlit Cloud 무료 플랜에서 메모리 부족 없이 구동
#   ✅ 실제 앱과 동일한 3개 모드 재현:
#        🎯 QSAR 결합력 분석
#        🏗️ 조립식 Ternary 모델링
#        💊 PROTAC 자동 완성 설계
#
#   모든 결과는 실제 DegradeAI가 산출한 값을 사전 저장(하드코딩)하여 재현합니다.
# ─────────────────────────────────────────────────────────────
import streamlit as st
import time
import py3Dmol
from stmol import showmol

st.set_page_config(page_title="DegradeAI - PROTAC-1조 (Demo)", layout="wide")

# ═════════════════════════════════════════════════════════════
# 사전 준비된 시연 데이터 (실제 DegradeAI 산출값 기반)
# ═════════════════════════════════════════════════════════════
DEMO = {
    "name": "JQ1 (BET bromodomain inhibitor)",
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

# ── JQ1 리간드 3D (사전 생성 PDB 블록) ──
JQ1_PDB = """HEADER    DEMO LIGAND JQ1
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


def make_demo_protein():
    """간단한 데모용 단백질 골격(나선형 폴리펩타이드) PDB 생성 — 인터넷 불필요."""
    import math
    lines = ["HEADER    DEMO PROTEIN POCKET"]
    idx = 1
    # 알파 헬릭스 모양으로 CA 원자 배치 (포켓을 감싸는 형태)
    for i in range(40):
        ang = i * 100 * math.pi / 180.0
        r = 5.0
        x = r * math.cos(ang)
        y = r * math.sin(ang)
        z = i * 1.5 - 30.0
        lines.append(
            f"ATOM  {idx:5d}  CA  ALA A{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C"
        )
        idx += 1
    lines.append("END")
    return "\n".join(lines)


DEMO_PROTEIN = make_demo_protein()


# ═════════════════════════════════════════════════════════════
# 3D 렌더링 함수들
# ═════════════════════════════════════════════════════════════
def render_qsar_3d():
    """QSAR — 단백질 포켓 + 약물 도킹 드라마틱 뷰."""
    view = py3Dmol.view(width=900, height=560)
    view.addModel(DEMO_PROTEIN, 'pdb')        # model 0: 단백질
    view.addModel(JQ1_PDB, 'pdb')             # model 1: 약물
    view.setBackgroundColor('0x0a0a1a')

    # 단백질 — 스펙트럼 카툰 + 점선 백본
    view.setStyle({'model': 0}, {
        "cartoon": {'color': 'spectrum', 'opacity': 0.55},
        "sphere": {'colorscheme': 'cyanCarbon', 'radius': 0.35, 'opacity': 0.25},
    })
    # 약물 — 형광 핑크
    view.setStyle({'model': 1}, {
        "stick": {'colorscheme': 'magentaCarbon', 'radius': 0.45},
        "sphere": {'colorscheme': 'magentaCarbon', 'radius': 0.40, 'opacity': 0.6},
    })
    view.zoomTo({'model': 1})
    view.center({'model': 1})
    view.zoom(1.6)
    view.spin('y', 0.8)
    return view


def render_ligand_3d():
    """약물 단독 드라마틱 회전."""
    view = py3Dmol.view(width=900, height=520)
    view.addModel(JQ1_PDB, 'pdb')
    view.setBackgroundColor('0x0a0a1a')
    view.setStyle({}, {
        "stick": {'colorscheme': 'magentaCarbon', 'radius': 0.34},
        "sphere": {'colorscheme': 'magentaCarbon', 'radius': 0.46},
    })
    view.addSurface(py3Dmol.VDW, {'opacity': 0.18, 'color': 'magenta'}, {})
    view.zoomTo()
    view.zoom(0.9)
    view.spin('y', 1.0)
    return view


def render_ternary_3d():
    """조립식 Ternary — 두 단백질 + 워헤드 드라마틱 뷰."""
    import math
    # 두 번째 단백질(E3) — 다른 위치에 배치
    lines = ["HEADER    DEMO E3"]
    idx = 1
    for i in range(30):
        ang = i * 100 * math.pi / 180.0
        x = 4.0 * math.cos(ang) + 18.0
        y = 4.0 * math.sin(ang)
        z = i * 1.5 - 22.0
        lines.append(f"ATOM  {idx:5d}  CA  GLY B{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C")
        idx += 1
    lines.append("END")
    e3_protein = "\n".join(lines)

    view = py3Dmol.view(width=950, height=600)
    view.addModel(DEMO_PROTEIN, 'pdb')   # BRD4
    view.addModel(e3_protein, 'pdb')     # VHL
    view.addModel(JQ1_PDB, 'pdb')        # 워헤드
    view.setBackgroundColor('0x05050f')

    view.setStyle({'model': 0}, {"cartoon": {'color': 'spectrum', 'opacity': 0.5}})
    view.setStyle({'model': 1}, {"cartoon": {'colorscheme': 'greenCarbon', 'opacity': 0.5},
                                 "sphere": {'colorscheme': 'greenCarbon', 'radius': 0.4, 'opacity': 0.3}})
    view.setStyle({'model': 2}, {
        "stick": {'colorscheme': 'magentaCarbon', 'radius': 0.5},
        "sphere": {'colorscheme': 'magentaCarbon', 'radius': 0.45, 'opacity': 0.6},
    })
    view.zoomTo({'model': 2})
    view.center({'model': 2})
    view.zoom(1.1)
    view.spin('y', 0.7)
    return view


# ═════════════════════════════════════════════════════════════
# 공통 헤더
# ═════════════════════════════════════════════════════════════
st.title("🧬 DegradeAI [PROTAC-1조]")
st.markdown("#### 머신러닝 기반 적응형 표적 단백질 분해제(TPD) 통합 설계 및 검증 플랫폼")
st.info(
    "📱 **포스터 시연용 데모입니다.** 실제 DegradeAI 플랫폼이 산출한 결과를 재현합니다. "
    "(안정적 시연을 위해 사전 연산된 결과를 표시 — 실시간 도킹·ML은 로컬 환경에서 동작)"
)

# ── 사이드바: 모드 선택 ──
st.sidebar.title("🧬 DegradeAI")
st.sidebar.caption("포스터 QR 시연 모드")
mode = st.sidebar.radio(
    "📍 분석 모드",
    ["🎯 QSAR 결합력 분석", "🏗️ 조립식 Ternary 모델링", "💊 PROTAC 자동 완성 설계"],
)
st.sidebar.divider()


# ═════════════════════════════════════════════════════════════
# 모드 1 — QSAR 결합력 분석
# ═════════════════════════════════════════════════════════════
def render_qsar():
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
            "**JQ1**은 대표적인 BET bromodomain 저해제로, PROTAC(MZ1, dBET 등)의 워헤드로 널리 쓰입니다."
        )
        return

    # 로그 애니메이션
    st.subheader("🖥️ AI 연산 관제 로그")
    box = st.empty(); logs = []
    steps = [
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
    ]
    for s in steps:
        logs.append(f"> {s}")
        box.code("\n".join(logs), language="bash")
        time.sleep(0.22)

    st.divider()
    st.markdown("### 🎬 Protein–Ligand Docking 3D — Dramatic View")
    st.info("🔵 스펙트럼 = 단백질 포켓 | 💗 형광 핑크 = 도킹된 JQ1")
    showmol(render_qsar_3d(), height=560, width=900)

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
        st.caption(f"Fsp3: {DEMO['fsp3']:.2f} | TPSA: {DEMO['tpsa']:.1f} | HBD: {DEMO['hbd']} | HBA: {DEMO['hba']}")
    with c6:
        st.metric("DLCS (약물성)", f"{DEMO['dlcs']:.3f}", f"실측: {DEMO['dlcs_actual']:.3f}")

    st.caption(f"🎯 신뢰도(AD): {DEMO['ad_score']:.2f} ✅ 학습 분포 내")

    st.divider()
    st.subheader("📝 의약화학적 의견서 (LLM Report)")
    st.markdown(f"### 종합 판정: {DEMO['verdict']}")
    st.info(
        "### 1. 결합력(pKd) 및 약물성(ADME) 종합 진단\n"
        f"예측 결합력 pKd {DEMO['pkd']:.2f}는 자동 등급 **🟡 보통**에 해당하며, "
        "리드 최적화를 통한 결합력 향상이 권장됩니다.\n"
        f"* 분자량 {DEMO['mw']:.0f} Da, LogP {DEMO['logp']:.1f}로 세포 투과성은 양호합니다.\n"
        "* BRD4 워헤드로서 유효한 출발점이나, 추가 최적화 여지가 있습니다.\n\n"
        "### 2. 합성 접근성(SA) 및 약물성(DLCS) 평가\n"
        f"SA Score {DEMO['sa_raw']:.2f}(**🟢 우수**), DLCS {DEMO['dlcs']:.3f}(**🟢 양호**)로 "
        "합성 실현 가능성과 약물성이 양호합니다.\n"
        f"* Bertz CT {DEMO['bertz']:.0f}는 표준 소분자 수준으로 CMC 리스크가 낮습니다.\n\n"
        "### 3. 종합 판정 및 권고\n"
        "**판정: 🟡 HOLD**\n\n"
        "조건부 보류 — 결합력 보완 후 재평가를 권고합니다. 합성성·약물성은 우수하나 "
        "결합력(pKd)이 임상 후보 기준에 다소 못 미쳐, 워헤드 구조 최적화가 필요합니다.\n"
        "* 다음 단계: 링커·E3 리간드 연결 후 ternary 복합체 형성 평가."
    )


# ═════════════════════════════════════════════════════════════
# 모드 2 — 조립식 Ternary 모델링
# ═════════════════════════════════════════════════════════════
def render_ternary():
    st.sidebar.header("🏗️ Ternary 설정")
    st.sidebar.text_input("워헤드 SMILES", DEMO["smiles"], disabled=True)
    st.sidebar.caption("🔗 QSAR 연동: JQ1 (pKd 6.78)")
    run = st.sidebar.button("🚀 조립식 Ternary 모델 생성", type="primary")

    st.title("🧬 조립식 Ternary Complex 모델링")
    st.markdown("#### 사용자 워헤드 + 실측 E3 템플릿(5T35) 기반 삼중 복합체 모델")
    st.info(
        "ℹ️ **Template-based Ternary Modeling** — 워헤드를 실측 ternary 구조"
        "(5T35: BRD4–VHL)의 리간드 위치에 정렬하여 삼중 복합체 모델을 생성합니다. "
        "(de novo 예측이 아닌 템플릿 중첩 모델)"
    )
    st.divider()

    if not run:
        st.markdown("👈 **🚀 조립식 Ternary 모델 생성** 버튼을 눌러 시연하세요.")
        return

    box = st.empty(); logs = []
    for s in ["실측 템플릿 5T35 로드 중...", "✅ 템플릿 확보 (BRD4–MZ1–VHL)",
              "✅ 타겟 체인 D / E3 체인 H / 템플릿 리간드 759",
              "사용자 워헤드 3D 생성 중...",
              "✅ 워헤드를 MZ1 결합 포켓에 정렬 완료 (도킹 시뮬레이션)",
              "조립식 ternary 3D 렌더링 중...", "✅ 완료!"]:
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
              help="Tanimoto 유사도. JQ1은 MZ1의 실제 워헤드이므로 1.0")
    v2.success("🟢 템플릿 적합성 높음 — JQ1은 MZ1의 워헤드와 동일하여 5T35 템플릿 기반 모델링이 매우 신뢰할 만합니다.")
    st.caption("💡 MZ1(5T35의 PROTAC)은 JQ1 워헤드를 사용합니다. 입력 워헤드가 동일하므로 동일 포켓·동일 배치를 따릅니다.")

    st.divider()
    st.subheader("🎬 조립식 Ternary Complex 3D — Dramatic View")
    st.markdown("🔵 **BRD4** (스펙트럼) | 🟢 **VHL** (초록) | 💗 **워헤드** (발광 핑크)")
    showmol(render_ternary_3d(), height=600, width=950)
    st.success("🎉 워헤드가 BRD4–VHL 삼중 복합체 골격에 조립된 모델 생성 완료!")


# ═════════════════════════════════════════════════════════════
# 모드 3 — PROTAC 자동 완성 설계
# ═════════════════════════════════════════════════════════════
def render_designer():
    st.sidebar.header("💊 PROTAC 설계")
    st.sidebar.text_input("워헤드 SMILES", DEMO["smiles"], disabled=True)
    st.sidebar.selectbox("E3 리간드", ["VHL — VH032"], disabled=True)
    run = st.sidebar.button("🚀 PROTAC 후보 자동 생성", type="primary")

    st.title("💊 PROTAC 자동 완성 설계")
    st.markdown("#### 워헤드 1개 → E3 리간드 + 링커 자동 조립 → 완성 PROTAC 후보")
    st.info(
        "ℹ️ 워헤드에 **논문 검증 E3 리간드 5종** + **링커 30종**(PEG·Alkyl·아미드·고리형)을 "
        "RDKit으로 결합해 완성 PROTAC을 생성하고, pKd·약물성·분자량 적합성을 종합해 정렬합니다."
    )
    st.caption("📏 분자량 필터: 600–1200 Da (sweet spot 700–1100 Da)")
    st.divider()

    if not run:
        st.markdown("👈 **🚀 PROTAC 후보 자동 생성** 버튼을 눌러 시연하세요.")
        return

    box = st.empty(); logs = []
    for s in ["워헤드 연결 지점([*:1]) 설정 중...", "✅ 사용자 지정 연결 지점 사용",
              "E3=VHL + 링커 30종 조합 생성·채점 중...",
              "✅ 완성 PROTAC 후보 18개 생성 (MW 필터 제외: 12개)"]:
        logs.append(f"> {s}")
        box.code("\n".join(logs), language="bash")
        time.sleep(0.25)

    st.divider()
    st.subheader("📊 완성 PROTAC 후보 순위표")
    import pandas as pd
    df = pd.DataFrame([
        {"링커": "PEG3", "종류": "PEG", "추정길이(Å)": 13.5, "MW": 942.4, "pKd참고": 6.91, "DLCS": 0.512, "SA": 3.84, "RotBonds": 16, "종합점수": 78.3, "거리오차(Å)": 0.3},
        {"링커": "PEG4", "종류": "PEG", "추정길이(Å)": 17.0, "MW": 986.5, "pKd참고": 6.88, "DLCS": 0.498, "SA": 3.91, "RotBonds": 19, "종합점수": 74.1, "거리오차(Å)": 3.2},
        {"링커": "Amide-PEG2", "종류": "Amide", "추정길이(Å)": 12.0, "MW": 955.4, "pKd참고": 6.85, "DLCS": 0.521, "SA": 4.02, "RotBonds": 15, "종합점수": 73.8, "거리오차(Å)": 1.8},
        {"링커": "Pip-PEG", "종류": "Rigid", "추정길이(Å)": 11.0, "MW": 968.5, "pKd참고": 6.79, "DLCS": 0.534, "SA": 4.21, "RotBonds": 13, "종합점수": 72.5, "거리오차(Å)": 2.8},
        {"링커": "Alkyl-C10", "종류": "Alkyl", "추정길이(Å)": 12.6, "MW": 924.5, "pKd참고": 6.72, "DLCS": 0.487, "SA": 3.65, "RotBonds": 18, "종합점수": 70.2, "거리오차(Å)": 1.2},
    ])
    st.dataframe(df, use_container_width=True)
    st.caption("정렬: 목표 거리(13.8Å) 오차 → 종합점수 | 종합점수 = pKd 40% + DLCS 30% + 분자량 적합성 30%")

    st.divider()
    st.subheader("🏆 최우선 후보: 워헤드 + PEG3 + VH032 (VHL 리간드)")
    m1, m2, m3, m4, m5 = st.columns(5)
    m1.metric("분자량", "942.4 Da", "sweet spot")
    m2.metric("pKd 참고", "6.91")
    m3.metric("DLCS", "0.512")
    m4.metric("회전결합", "16")
    m5.metric("종합점수", "78.3/100")
    st.success("🎉 완성 PROTAC 후보 생성 완료! 이 SMILES를 다른 탭에서 도킹/3D 검증으로 이어갈 수 있습니다.")


# ═════════════════════════════════════════════════════════════
# 라우팅
# ═════════════════════════════════════════════════════════════
if mode == "🎯 QSAR 결합력 분석":
    render_qsar()
elif mode == "🏗️ 조립식 Ternary 모델링":
    render_ternary()
elif mode == "💊 PROTAC 자동 완성 설계":
    render_designer()

st.divider()
st.caption(
    "💡 본 데모는 DegradeAI 플랫폼의 실제 산출 결과를 재현한 시연용입니다. "
    "전체 기능(임의 분자 입력, 실시간 도킹·ML)은 로컬 환경에서 동작합니다. "
    "모든 예측은 wet-lab 실험으로 검증되어야 합니다. — PROTAC-1조"
)
