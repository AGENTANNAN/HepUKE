# ============================================================
# Paper: Search for e+e- → χcJ π+π- and a charmonium-like
#        structure decaying to χcJ π± between 4.18 and 4.60 GeV
# arXiv: 2012.02682v2
# ============================================================
# χcJ (J=0,1,2) reconstructed via χcJ → γ J/ψ, J/ψ → ℓ+ℓ- (ℓ=e,μ).
# Same event selection for all three χcJ states; separation by
# photon energy in the π+π- recoil system performed in ROOT.

###
### Datasets: 15 energy points, BOSS 703 (Table I)
###
data_4180 = DatasetManager.real_data.find("703_4180")  # 4.178 GeV, 3194.0 pb-1
data_4190 = DatasetManager.real_data.find("703_4190")  # 4.189 GeV, 526.7 pb-1
data_4200 = DatasetManager.real_data.find("703_4200")  # 4.200 GeV, 526.0 pb-1
data_4210 = DatasetManager.real_data.find("703_4210")  # 4.210 GeV, 517.1 pb-1
data_4220 = DatasetManager.real_data.find("703_4220")  # 4.219 GeV, 514.6 pb-1
data_4230 = DatasetManager.real_data.find("703_4230")  # 4.226 GeV, 1056.4 pb-1
data_4240 = DatasetManager.real_data.find("703_4240")  # 4.236 + 4.244 GeV, 530.3 + 538.1 pb-1
data_4260 = DatasetManager.real_data.find("703_4260")  # 4.258 GeV, 828.4 pb-1
data_4270 = DatasetManager.real_data.find("703_4270")  # 4.267 GeV, 531.1 pb-1
data_4280 = DatasetManager.real_data.find("703_4280")  # 4.278 GeV, 175.7 pb-1
data_4360 = DatasetManager.real_data.find("703_4360")  # 4.358 GeV, 543.9 pb-1
data_4420 = DatasetManager.real_data.find("703_4420")  # 4.416 GeV, 1044.0 pb-1
data_4530 = DatasetManager.real_data.find("703_4530")  # 4.527 GeV, 112.1 pb-1
data_4600 = DatasetManager.real_data.find("703_4600")  # 4.600 GeV, 586.9 pb-1

data_points = [data_4180, data_4190, data_4200, data_4210, data_4220,
               data_4230, data_4240, data_4260, data_4270, data_4280,
               data_4360, data_4420, data_4530, data_4600]
inc_mc_points = data_points.map { |d| DatasetManager.inclusive_mc.find("#{d.boss}_#{d.sample_name}") }

# ============================================================
# Decay cards: one per χcJ state (J=0,1,2)
# χcJ → γ J/ψ, J/ψ → e+e- (μ+μ- equivalent; both combined in analysis)
# ============================================================
decay_card_c0 = <<~DECAYCARD
    Decay psi(4260)
    1.000 chi_c0 pi+ pi-  PHSP;
    Enddecay

    Decay chi_c0
    1.000 gamma J/psi  VSP_PWAVE;
    Enddecay

    Decay J/psi
    1.000 e+ e-  VLL;
    Enddecay

    End
DECAYCARD

decay_card_c1 = <<~DECAYCARD
    Decay psi(4260)
    1.000 chi_c1 pi+ pi-  PHSP;
    Enddecay

    Decay chi_c1
    1.000 gamma J/psi  VSP_PWAVE;
    Enddecay

    Decay J/psi
    1.000 e+ e-  VLL;
    Enddecay

    End
DECAYCARD

decay_card_c2 = <<~DECAYCARD
    Decay psi(4260)
    1.000 chi_c2 pi+ pi-  PHSP;
    Enddecay

    Decay chi_c2
    1.000 gamma J/psi  VSP_PWAVE;
    Enddecay

    Decay J/psi
    1.000 e+ e-  VLL;
    Enddecay

    End
DECAYCARD

# Exclusive MC for each χcJ state (multi-energy scan)
exMC_c0 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_chic0_pipi_ee"
  config.events        = 500_000
  config.decay_card    = decay_card_c0
  config.cross_section = :default
end

exMC_c1 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_chic1_pipi_ee"
  config.events        = 500_000
  config.decay_card    = decay_card_c1
  config.cross_section = :default
end

exMC_c2 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_chic2_pipi_ee"
  config.events        = 500_000
  config.decay_card    = decay_card_c2
  config.cross_section = :default
end

# ============================================================
# Algorithm: χcJ π+π- (single selection for all J=0,1,2)
# χcJ separation via γ energy in π+π- recoil frame (ROOT level)
# ============================================================
alg = Algorithm.new("ChicJPipPim")
alg.set_header(["ChicJPipPimAlg/ChicJPipPim.h"])
alg.set_constant({})

sel = Selection.new

# --- Charged track selection ---
# 4 tracks: π+π- + ℓ+ℓ-, net charge zero
sel.select_track do
  cos_theta   0.93
  Vz          10.0
  Vr          1.0
  nChrp       "==2"
  nChrn       "==2"
  nNet        "==0"
end

# --- Photon selection ---
# At least 1 photon (from χcJ → γ J/ψ)
# Barrel: E > 25 MeV, |cosθ|<0.80; Endcap: E > 50 MeV, 0.86<|cosθ|<0.92
# EMC time within 700 ns; > 20° from nearest charged track
sel.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    20.0
  nGam              ">=1"
end

# --- Particle identification ---
# Pion candidates: tracks with p < 1.0 GeV/c
# Lepton candidates: tracks with p > 1.0 GeV/c
# Within leptons: E_EMC/P_MDC > 0.7c → electron; < 0.3c → muon
# Both e+e- and μ+μ- final states are combined in the analysis
sel.pid(method: :probability) do
  prob_cut 0.001
  identify :electron, against: [:pion, :kaon, :proton]
  nep "==1"
  nem "==1"
end

sel.assign({chrgp: :pip, chrgn: :pim})

# --- 5C kinematic fit: 4C + J/ψ mass constraint ---
# J/ψ → e+e- nominal hypothesis; μ+μ- hypothesis handled separately
sel.kinematic_fit([:pip, :pim, :ep, :em, :gamma]) do
  nominal
  constrain_four_momentum
  invariant_mass_of(:ep, :em).constrain_to_nominal_mass_of(:"J/psi")
  chi2_cut 50
end

# --- Competing μ+μ- hypothesis ---
sel.assign({ep: :mup, em: :mum})
  .kinematic_fit([:pip, :pim, :mup, :mum, :gamma]) do
    use_track_index_from_nominal_kmfit
    constrain_four_momentum
    invariant_mass_of(:mup, :mum).constrain_to_nominal_mass_of(:"J/psi")
  end

alg.note(:chicJ_separation, "χcJ (J=0,1,2) separated by photon energy in π+π- recoil rest frame (ROOT level); Fig. 2 shows selection windows")
alg.note(:lepton_id_momentum, "Pion/lepton separation via p<1.0 GeV/c (pion) vs p>1.0 GeV/c (lepton); e/μ separation via E_EMC/P_MDC > 0.7c (e) vs < 0.3c (μ)")
alg.note(:best_candidate, "Best candidate selected by lowest χ²_5C; multiple candidates rare in practice")
alg.note(:both_lepton_channels, "Both J/ψ→e+e- and J/ψ→μ+μ- combined; combined branching fraction B(J/ψ→ℓ+ℓ-) used in cross-section calculation")
alg.note(:background_vetoes, "ROOT-level vetoes: cos(α_π+π-) < 0.98 (Bhabha), m_rec veto for η/η' J/ψ, ω veto via χcJ recoil mass, ψ(2S)/X(3872) veto via m(π+π-J/ψ)")
alg.note(:isr_correction, "ISR correction factor κ = 0.64 (conservative) applied in ROOT; vacuum polarization factor 1/|1-Π(s)|² from Ref. [41]")
alg.note(:no_signal, "No significant signal observed at any energy; 90% CL upper limits set for all three χcJ channels")
alg.note(:systematics, "Total systematic uncertainty 4.7-11.0% depending on χcJ mode and energy; Gaussian smearing of efficiency in UL calculation")

alg.with_decay_card(decay_card_c1)  # representative χcJ=1 decay card
   .apply(sel)

alg.execute_on(data_points + inc_mc_points + exMC_c0 + exMC_c1 + exMC_c2)