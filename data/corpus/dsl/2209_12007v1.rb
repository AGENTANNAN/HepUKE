# Paper: 2209.12007v1 — Search for X(3872) via e⁺e⁻ → π⁺π⁻J/ψ cross section
# Ordinary analysis at 4 energy points (3.808–3.896 GeV)
# J/ψ reconstructed via e⁺e⁻ and μ⁺μ⁻ channels

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_3810 = DatasetManager.real_data.find("703_3810")
data_3872 = DatasetManager.real_data.find("703_3872")
data_3900 = DatasetManager.real_data.find("703_3900")

all_data = [data_3810, data_3872, data_3900]

incMC_3810 = DatasetManager.inclusive_mc.find("703_3810")
incMC_3872 = DatasetManager.inclusive_mc.find("703_3872")
incMC_3900 = DatasetManager.inclusive_mc.find("703_3900")

all_incMC = [incMC_3810, incMC_3872, incMC_3900]

# Exclusive MC: e⁺e⁻ → π⁺π⁻J/ψ, J/ψ → e⁺e⁻ (representative; μμ channel similar)
decay_card = <<~DECAYCARD
Decay psi(4260)
1 pi+ pi- J/psi PHSP;
Enddecay
Decay J/psi
1 e+ e- PHSP;
Enddecay
End
DECAYCARD

sig_mc = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "sig_pipi_jpsi"
  config.events        = 200_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

alg = Algorithm.new("X3872PipiJpsi")
alg.set_header(["X3872PipiJpsiAlg/X3872PipiJpsi.h"])
    .set_constant({ "ECMS" => [:double, 3.808] })
    .with_decay_card(decay_card)

event_selection = Selection.new
  .select_track do
    cos_theta  0.93
    Vz         10.0
    Vr         1.0
    nChrp      "==2"
    nChrn      "==2"
    nNet       "==0"
  end
  .pid do
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 1.1
    nlp "==1"
    nlm "==1"
  end
  .remove([:lp <= :chrgp, :lm <= :chrgn])
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip "==1"
    npim "==1"
  end
  .kinematic_fit([:pip, :pim, :lp, :lm]) do
    invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
    constrain_four_momentum
    chi2_cut 60
    nominal
  end

alg.note(:energy_points, "Paper uses 4 energy points at 3807.7, 3867.4, 3871.3, 3896.2 MeV. Dataset table provides 703_3810 (3807.65), 703_3872 (3867.41), 703_3900 (3896.24). The 3871.3 MeV point is within the 703_3872 sample; BOSS treats it as the same dataset.")
alg.note(:lepton_channels, "Paper uses J/ψ→e⁺e⁻ and J/ψ→μ⁺μ⁻ channels. DSL uses identify_high_momentum_leptons which classifies e vs μ by EMC energy (electron: EMC>1.1 GeV, muon: EMC<0.35 GeV). The requirement that both leptons be same type and the channel-specific analysis are ROOT-level.")
alg.note(:muon_channel, "Muon channel (J/ψ→μ⁺μ⁻) requires separate MC sample and analysis. The muon EMC<0.35 GeV cut is not DSL-tunable in identify_high_momentum_leptons; this cut is applied at ROOT level.")
alg.note(:momentum_cuts, "Paper: low momentum (<0.6 GeV/c) → pion, high momentum (>1.0 GeV/c) → lepton. DSL identify_high_momentum_leptons handles the 1.0 GeV/c threshold; pion momentum cut is ROOT-level.")
alg.note(:photon_conversion_veto, "Paper: cosθ_ππ<0.95, cosθ_πe<0.98 to suppress photon conversion background. ROOT-level cut.")
alg.note(:track_region_gap, "Paper leaves 0.6–1.0 GeV/c gap unclassified. DSL handles this naturally: tracks < 1.0 GeV/c remain after lepton removal and become pions.")
alg.note(:x3872_search, "Paper searches for X(3872) in M(π⁺π⁻J/ψ) spectrum. ROOT analysis beyond DSL scope.")

alg.with_decay_card(decay_card).apply(event_selection)
alg.execute_on(all_data + all_incMC + sig_mc)