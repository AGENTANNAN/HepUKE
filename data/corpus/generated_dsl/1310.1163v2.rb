# ================================================================
# e+e- -> pi (D D*) at 4.260 GeV  --  partial reconstruction
#   Mode 1 : e+e- -> pi+ D0  D*- , D0 -> K- pi+
#   Mode 2 : e+e- -> pi- D+  D*0 , D+ -> K- pi+ pi+
# Only the bachelor pion and one D meson are detected; the D* is
# inferred from energy-momentum conservation (recoil of pi + D).
# ================================================================

### ---------------------------------------------------------------
### Datasets : 4.260 GeV real data + matching inclusive MC
### ---------------------------------------------------------------
data_4260  = DatasetManager.real_data.find("703_4260")      # 4.260 GeV real data
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")  # corresponding inclusive MC

### ---------------------------------------------------------------
### Decay cards (EvtGen)
### ---------------------------------------------------------------
decay_card_mode1 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 pi+ D0 D*- PHSP;
    Enddecay

    Decay D0
    1.0000 K- pi+ PHSP;
    Enddecay

    Decay D*-
    1.0000 anti-D0 pi- PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_mode2 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 pi- D+ D*0 PHSP;
    Enddecay

    Decay D+
    1.0000 K- pi+ pi+ PHSP;
    Enddecay

    Decay D*0
    1.0000 D0 pi0 PHSP;
    Enddecay

    Decay D0
    1.0000 K- pi+ PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### ---------------------------------------------------------------
### Exclusive MC : 100k events per mode
### ---------------------------------------------------------------
exMC_mode1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4260_piplus_D0_Dstar"
  config.related_dataset = data_4260
  config.events          = 100_000
  config.decay_card      = decay_card_mode1
  config.cross_section   = :default
end

exMC_mode2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4260_piminus_Dplus_Dstar0"
  config.related_dataset = data_4260
  config.events          = 100_000
  config.decay_card      = decay_card_mode2
  config.cross_section   = :default
end

### ---------------------------------------------------------------
### BOSS event selection -- Mode 1 : pi+ D0 D*- , D0 -> K- pi+
### ---------------------------------------------------------------
alg_mode1 = Algorithm.new("PiD0Dstar")
alg_mode1.set_header(["PiD0DstarAlg/PiD0Dstar.h"])
         .set_constant({ "ECMS" => [:double, 4.260] })
         .set_alias({ "std::vector<double>" => "Vdouble" })

sel_mode1 = Selection.new
sel_mode1
  .select_track {
      cos_theta 0.93        # |cos(theta)| < 0.93
      Vz        10.0        # |Vz| < 10 cm
      Vr        1.0         # Vr < 1 cm
      nTot      ">=3"       # at least three good charged tracks
  }
  .pid(method: :probability) {
      prob_cut 0.001                              # PID probability > 0.001
      identify :kaon, against: [:pion, :proton]   # K+/- separated from pi and p
      identify :pion, against: [:kaon, :proton]   # pi+/- separated from K and p
      nkm  ">=1"                                  # at least one K-
      npip ">=2"                                  # at least two pi+
  }
  # Mass-constrained reconstruction of D0 -> K- pi+ (loose chi2 default;
  # published tight chi2<30 for the 2-constraint fit is applied later in ROOT)
  .kalman_kinematic_fit([:km, :pip]) {
      invariant_mass_of(:km, :pip).within(1.8498, 1.8798)          # |M(K-pi+) - m(D0)| < 15 MeV/c^2
      invariant_mass_of(:km, :pip).constrain_to_nominal_mass_of(:D0)
      chi2_cut 200
      nD0 ">=1"
  }
  # Partial reconstruction: D*- is NOT built from its decay products but inferred
  # from energy-momentum conservation. In decay_card_mode1 recID 3 is D*- ; its
  # daughters are expanded automatically. D0 (recID 2) is taken from the fit above.
  .partial_miss([3])

alg_mode1.note(:background_veto, "post-fit M(pi+ D0) > 2.02 GeV/c^2 veto on the pi+ D0 tag to reject the D*+ D*- background; applied on the ROOT side after this BOSS selection")

alg_mode1.with_decay_card(decay_card_mode1).apply(sel_mode1)
alg_mode1.execute_on([data_4260, incMC_4260, exMC_mode1])

### ---------------------------------------------------------------
### BOSS event selection -- Mode 2 : pi- D+ D*0 , D+ -> K- pi+ pi+
### ---------------------------------------------------------------
alg_mode2 = Algorithm.new("PiDpDstar0")
alg_mode2.set_header(["PiDpDstar0Alg/PiDpDstar0.h"])
         .set_constant({ "ECMS" => [:double, 4.260] })
         .set_alias({ "std::vector<double>" => "Vdouble" })

sel_mode2 = Selection.new
sel_mode2
  .select_track {
      cos_theta 0.93        # |cos(theta)| < 0.93
      Vz        10.0        # |Vz| < 10 cm
      Vr        1.0         # Vr < 1 cm
      nTot      ">=3"       # at least three good charged tracks
  }
  .pid(method: :probability) {
      prob_cut 0.001                              # PID probability > 0.001
      identify :kaon, against: [:pion, :proton]   # K+/- separated from pi and p
      identify :pion, against: [:kaon, :proton]   # pi+/- separated from K and p
      nkm  ">=1"                                  # at least one K-
      npip ">=2"                                  # at least two pi+
      npim ">=1"                                  # pi- D+ tag requires at least one pi-
  }
  # Mass-constrained reconstruction of D+ -> K- pi+ pi+
  .kalman_kinematic_fit([:km, :pip, :pip]) {
      invariant_mass_of(:km, :pip, :pip).within(1.8547, 1.8847)    # |M(K-pi+pi+) - m(D+)| < 15 MeV/c^2
      invariant_mass_of(:km, :pip, :pip).constrain_to_nominal_mass_of(:Dplus)
      chi2_cut 200
      nDplus ">=1"
  }
  # Partial reconstruction: D*0 inferred from energy-momentum conservation.
  # In decay_card_mode2 recID 3 is D*0 ; its daughters are expanded automatically.
  .partial_miss([3])

alg_mode2.with_decay_card(decay_card_mode2).apply(sel_mode2)
alg_mode2.execute_on([data_4260, incMC_4260, exMC_mode2])