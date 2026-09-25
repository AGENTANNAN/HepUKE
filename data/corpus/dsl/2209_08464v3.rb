# Paper: 2209.08464v3 — PWA of Λc⁺ → Λπ⁺π⁰
# Tag-based single-tag analysis at 7 energy points (4.600–4.699 GeV)

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_4600 = DatasetManager.real_data.find("703_4600")
data_4610 = DatasetManager.real_data.find("706_4610")
data_4620 = DatasetManager.real_data.find("706_4620")
data_4640 = DatasetManager.real_data.find("706_4640")
data_4660 = DatasetManager.real_data.find("706_4660")
data_4680 = DatasetManager.real_data.find("706_4680")
data_4700 = DatasetManager.real_data.find("706_4700")

all_data = [data_4600, data_4610, data_4620, data_4640, data_4660, data_4680, data_4700]

incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_4610 = DatasetManager.inclusive_mc.find("706_4610")
incMC_4620 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4640 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")

all_incMC = [incMC_4600, incMC_4610, incMC_4620, incMC_4640, incMC_4660, incMC_4680, incMC_4700]

# Exclusive signal MC: e⁺e⁻ → Λc⁺Λc⁻ with Λc⁺ → Λπ⁺π⁰
decay_card = <<~DECAYCARD
Decay psi(4260)
1 Lambda_c+ anti-Lambda_c- PHSP;
Enddecay
Decay Lambda_c+
1 Lambda pi+ pi0 PHSP;
Enddecay
Decay anti-Lambda_c-
1 anti-Lambda pi- pi0 PHSP;
Enddecay
Decay Lambda
1 p+ pi- PHSP;
Enddecay
Decay anti-Lambda
1 anti-p- pi+ PHSP;
Enddecay
End
DECAYCARD

sig_mc = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "sig_lambdac_lambda_pi_pi0"
  config.events        = 1_000_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

alg = TagAnalysis.new("LambdacToLambdaPiPi0")
alg.set_header(["LambdacToLambdaPiPi0Alg/LambdacToLambdaPiPi0.h"])
    .set_constant({ "ECMS" => [:double, 4.600] })
    .with_decay_card(decay_card)

# Tag side: anti-Λc⁻ tagged via hadronic decays
# Paper uses 10 hadronic decay modes; here we declare the most common ones.
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKPiP
  t.charm -1
end

# Signal side: Λc⁺ → Λπ⁺π⁰ (Λ → pπ⁻, π⁰ → γγ)
alg.signal_side do |s|
  s.photons 2
  s.charged(prp: 1, pim: 1, pip: 1)
  s.require_charge 1
end

# Kinematic fit: 4-momentum conservation + Λ mass + π⁰ mass + tag Λc⁺ mass
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:prp, :pim).constrain_to_nominal_mass_of(:Lambda)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

# Inexpressible BOSS procedures
alg.note(:tag_mode_list, "Paper uses 10 hadronic Λc decay modes for ST; only :LambdacPtoKPiP declared. Other 9 modes require verification against tag_modes.yaml and may need additional mode entries.")
alg.note(:track_quality_cuts, "Paper applies |cosθ|<0.93, Vz<10cm (20cm for Λ daughters), Vr<1cm for charged tracks. These are BOSS defaults except Vz for Λ daughters — handled by TagAnalysis framework.")
alg.note(:photon_cuts, "Paper: photon E>25MeV barrel, E>50MeV endcap, EMC time within 700ns. π⁰: M(γγ)∈[0.115,0.150]GeV with 1C Kalman fit constraining to π⁰ mass. Λ: vertex fit χ²<100, decay length>2σ, M(pπ⁻)∈[1.111,1.121]GeV.")
alg.note(:signal_side_deltaE, "Paper: Λc⁺ candidate with min|ΔE| retained, ΔE∈[-0.03,0.02]GeV. Applied at analysis level; tag-side mBC and ΔE stored unconditionally per store-not-cut rule.")
alg.note(:sigma0_veto, "Paper: Σ⁰ veto — reject events satisfying both Λc⁺→Λπ⁺π⁰ and Λc⁺→Σ⁰π⁺ selection. Applied at ROOT analysis level.")
alg.note(:three_c_fit, "Paper describes a 3C kinematic fit (π⁰ mass, Λ mass, recoil mass). DSL emits a 7C fit (4C + 3 mass constraints) with equivalent physics content.")
alg.note(:pwa, "Partial wave analysis (PWA) of Λc⁺→Λπ⁺π⁰ including ρ(770)⁺, Σ(1385)⁺, Σ(1385)⁰ components is performed in ROOT analysis beyond the DSL scope.")

alg.apply
alg.execute_on(all_data + all_incMC + sig_mc)