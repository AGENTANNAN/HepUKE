# BOSS Ruby DSL for BESIII paper 1802.04057v2
# Search for ψ(3686) → Λc+ pbar e+e- (FCNC decay)
# Data: ψ(3686), 448.1×10⁶ events
# Reconstructed via Λc+ → p K- π+
# Final state: p pbar K- π+ e+ e- (6 charged tracks, zero net charge)
# Vertex fit + 4C kinematic fit, χ² < 200
# No signal observed; upper limit on BF at 90% CL: 1.7×10⁻⁶

### Dataset ###
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card: ψ(3686) → Λc+ pbar e+ e-, Λc+ → p K- π+
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.000  Lambda_c+  anti-p-  e+  e-  PHSP;
    Enddecay

    Decay Lambda_c+
    1.000  p+  K-  pi+  PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC sample
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_psip_Lambdac_pbar_ee"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card_signal
  config.cross_section = :default
end

### Event selection ###
# 6 charged tracks: p, pbar, K-, π+, e+, e-
# |cosθ| < 0.93, |Vz| < 10 cm, Vr < 1 cm, net charge = 0
# PID: highest confidence level assigns type. No extra tracks allowed (nTot "==6").
# Vertex fit must converge.
# 4C kinematic fit: χ² < 200.
# Anti-Lambda veto: M(pbar π+) > 1.13 GeV/c².

alg = Algorithm.new("PsipLambdacPbarEE")
alg.set_header(["PsipLambdacPbarEEAlg/PsipLambdacPbarEE.h"])
   .set_constant({"ECMS" => [:double, 3.686]})

sel = Selection.new
sel.select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nTot "==6"
    nNet "==0"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    nprp ">=2"
    nkm ">=1"
    npip ">=1"
    nlp ">=1"; nlm ">=1"
  end
  .kinematic_fit([:prp, :prm, :km, :pip, :lp, :lm]) do
    nominal
    vertex_fit([0, 1, 2, 3, 4, 5])
    constrain_four_momentum
    chi2_cut 200
  end

alg
  .note(:pid_correction_method, "track reconstruction and PID: 1.0% per track; lepton tracking/PID: 2.5% per lepton (low momentum)")
  .note(:helix_correction, "4C kinematic fit efficiency difference 1.0% from control sample ψ(3686)→π+π-J/ψ→π+π-ppbarπ+π-")
  .note(:background_veto, "anti-Lambda veto: M(pbar π+) > 1.13 GeV/c² to reject ψ(3686)→γχcJ, χcJ→pK-antiLambda backgrounds with γ conversion")
  .note(:efficiency_curve, "Λc+ signal region M(pK-π+) in [2.25,2.32] GeV/c²; sidebands [2.06,2.23]+[2.34,2.40] GeV/c²; VMD signal model with 34.3% model uncertainty")
  .with_decay_card(decay_card_signal)
  .apply(sel)

### Execute ###
alg.execute_on([psip_data, psip_incMC, exMC_signal])