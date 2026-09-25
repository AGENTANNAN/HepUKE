# ============================================================
# Paper: 2005.05850v1
# Measurement of Born Cross Sections for:
#   (1) e+e- -> Ds+ Ds1(2460)- + c.c.
#   (2) e+e- -> Ds*+ Ds1(2460)- + c.c.
# ============================================================

### Common datasets ###
data_4470 = DatasetManager.real_data.find("703_4470")   # sqrt(s) ~ 4.467 GeV
data_4530 = DatasetManager.real_data.find("703_4530")   # sqrt(s) ~ 4.527 GeV
data_4575 = DatasetManager.real_data.find("703_4575")   # sqrt(s) ~ 4.575 GeV
data_4600 = DatasetManager.real_data.find("703_4600")   # sqrt(s) ~ 4.600 GeV

incMC_4470 = DatasetManager.inclusive_mc.find("703_4470")
incMC_4530 = DatasetManager.inclusive_mc.find("703_4530")
incMC_4575 = DatasetManager.inclusive_mc.find("703_4575")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")

# NOTE: Energy points 4.550, 4.560, 4.570, 4.580, 4.590 GeV are not in the
# BES3_dataset table. They were small-statistics scan points (~8 pb^-1 each)
# used only for upper-limit determination. Excluded from execute_on.

# ============================================================
# Process 1: e+e- -> Ds+ Ds1(2460)- + c.c.
# Reconstruct Ds+ -> K+K-pi+, identify Ds1(2460)- via recoil mass
# ============================================================

decay_card_ds_ds1 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D_s+ D_s1-  PHSP;
    Enddecay
    Decay D_s+
    1.0000 K+ K- pi+   PHSP;
    Enddecay
    Decay D_s1-
    1.0000 D_s*+ pi0   PHSP;
    Enddecay
    Decay D_s*+
    1.0000 gamma D_s+   VSP_PWAVE;
    Enddecay
    Decay pi0
    1.0000 gamma gamma  PHSP;
    Enddecay
    End
DECAYCARD

exMC_ds_ds1 = DatasetManager.create_exclusive_mc_for([data_4470, data_4530, data_4575, data_4600]) do |config|
  config.sample_name   = "sig_DsDs1_2460"
  config.events        = 200000
  config.decay_card    = decay_card_ds_ds1
  config.cross_section = :default
end

alg_ds_ds1 = Algorithm.new("DsDs1_2460")
alg_ds_ds1.set_header(["DsDs1_2460Alg/DsDs1_2460.h"])

sel_ds_ds1 = Selection.new
sel_ds_ds1.select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"       # K+ and pi+
    nChrn ">=1"       # K-
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp ">=1"
    nkm ">=1"
    npip ">=1"
  end
  # Partial reconstruction: reconstruct D_s+ (recID 1), Ds1 inferred from recoil
  # Ds1(2460)- mass = 2.4595 GeV/c^2; identified via ROOT-level fit to M_rec+ spectrum
  .partial_rec([1]) do
    best_combination_by_mass :"D_s+", 1.96835
  end

alg_ds_ds1
  .note(:tag_mode_unavailable, "Ds+ reconstructed via two sub-modes:
    (1) Ds+ -> phi(->K+K-)pi+ with |M(K+K-)-m_phi|<15 MeV/c^2;
    (2) Ds+ -> K*(892)0(->K-pi+)K+ with |M(K-pi+)-m_K*0|<84 MeV/c^2.
    Both sub-modes share the K+K-pi+ final state; combinatorial reconstruction
    uses all three tracks. Invariant-mass windows on intermediate resonances
    cannot be expressed inside partial_rec DSL block.")
  .note(:recoil_mass_fitting, "Ds1(2460)- identified via recoil mass
    M_rec(Ds+) = M_recoil(KKpi) + M(KKpi) - m_Ds+. Unbinned ML fit to
    M(KKpi) in 4.0 MeV/c^2 M_rec+ bins extracts Ds+ signal yields; Ds1
    yield extracted from fit to M_rec+ distribution with MC-derived signal
    shape and polynomial background.")
  .note(:cross_section_measurement, "Born cross section:
    sigma_B = Nfit / (Lint * (1+delta) * (1+delta_vp) * eps_Ds)
    where eps_Ds = eps * B(Ds+->K+K-pi+). ISR factor (1+delta) from KKMC
    QED calculation with 1% accuracy. Vacuum polarization (1+delta_vp) near
    0.055 from Ref.[41]. Upper limits at 90% C.L. when significance < 3sigma.")
  .note(:systematic_uncertainties, "Multiplicative: tracking/PID 3.5%, MC stats
    0.5%, ISR correction 4.6-13.1%, luminosity 0.7-0.8%, BFs 3.2%, MC generator
    1.3%. Additive from fit: Ds+ mass resolution, bin width, Ds1 mass,
    background shape, fit range. Total 8.3-22.3% at significant points.")
  .with_decay_card(decay_card_ds_ds1)
  .apply(sel_ds_ds1)

# ============================================================
# Process 2: e+e- -> Ds*+ Ds1(2460)- + c.c.
# Reconstruct Ds*+ -> gamma Ds+, Ds+ -> K+K-pi+, Ds1 via recoil
# ============================================================

decay_card_dsstar_ds1 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D_s*+ D_s1-  PHSP;
    Enddecay
    Decay D_s*+
    1.0000 gamma D_s+    VSP_PWAVE;
    Enddecay
    Decay D_s+
    1.0000 K+ K- pi+    PHSP;
    Enddecay
    Decay D_s1-
    1.0000 D_s*+ pi0    PHSP;
    Enddecay
    Decay pi0
    1.0000 gamma gamma   PHSP;
    Enddecay
    End
DECAYCARD

exMC_dsstar_ds1 = DatasetManager.create_exclusive_mc_for([data_4600]) do |config|
  config.sample_name   = "sig_DsStarDs1_2460"
  config.events        = 100000
  config.decay_card    = decay_card_dsstar_ds1
  config.cross_section = :default
end

alg_dsstar_ds1 = Algorithm.new("DsStarDs1_2460")
alg_dsstar_ds1.set_header(["DsStarDs1_2460Alg/DsStarDs1_2460.h"])

sel_dsstar_ds1 = Selection.new
sel_dsstar_ds1.select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"
    nChrn ">=1"
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 20.0
    nGam ">=1"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp ">=1"
    nkm ">=1"
    npip ">=1"
  end
  # Partial reconstruction: reconstruct D_s*+ (recID 1) -> gamma + D_s+ (->K+K-pi+)
  # Ds1(2460)- inferred from recoil against D_s*+
  .partial_rec([1]) do
    best_combination_by_mass :"D_s*+", 2.1122
  end

alg_dsstar_ds1
  .note(:two_c_kinematic_fit, "2C mass-constrained fit to nominal Ds+ and Ds*+
    masses applied before recoil mass study. chi2_2C < 10 required. The
    2C fit improves the D_s*+ recoil mass resolution; cannot be expressed
    inside partial_rec block as partial_rec replaces kinematic_fit entirely.")
  .note(:tag_mode_unavailable, "Ds+ sub-modes (phi-pi+ and K*0-K+) share the
    same K+K-pi+ final state. Intermediate resonance mass windows applied
    before 2C fit: |M(K+K-)-m_phi|<15 MeV/c^2 or |M(K-pi+)-m_K*0|<84 MeV/c^2.")
  .note(:recoil_mass_fitting, "Ds1(2460)- signal extracted from unbinned ML fit
    to M_rec(Ds*+) distribution; signal modeled by Crystal Ball function,
    background by ARGUS function. Upper limit at 90% C.L. for sqrt(s)=4.590 GeV.")
  .note(:cross_section_measurement, "Born cross section at sqrt(s)=4.600 GeV.
    Significance 5.9 sigma at 4.600 GeV, 2.0 sigma at 4.590 GeV.
    Correction factors (1+delta), (1+delta_vp) from KKMC and Ref.[41].")
  .with_decay_card(decay_card_dsstar_ds1)
  .apply(sel_dsstar_ds1)

# ============================================================
# Execute
# ============================================================

alg_ds_ds1.execute_on([data_4470, data_4530, data_4575, data_4600,
                       incMC_4470, incMC_4530, incMC_4575, incMC_4600] + exMC_ds_ds1)

alg_dsstar_ds1.execute_on([data_4600, incMC_4600] + exMC_dsstar_ds1)