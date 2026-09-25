DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_4600  = DatasetManager.real_data.find("703_4600")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")

decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.000 p+ K_S0 eta PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.000 anti-p- pi+ PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Lc_pKs0eta"
  config.related_dataset = data_4600
  config.events          = 1_000_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

alg = Algorithm.new("LcPpKsEta")
alg.set_header(["LcPpKsEtaAlg/LcPpKsEta.h"])
   .set_constant({ "ECMS" => [:double, 4.59953] })

sel = Selection.new
sel.select_track {
      nChrp ">=1"
      nChrn ">=3"
      nNet "==0"
      cos_theta 0.93
      Vz 10.0
      Vr 1.0
    }
   .select_photon {
      tdc_emc_start 0
      tdc_emc_end 700
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam ">=2"
    }
   .pid(method: :probability) {
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]
      nprp ">=1"
    }
   .remove([:prp <= :chrgp])
   .assign({:chrgp => :pip, :chrgn => :pim})
   .secondary_vertex_fit([:pip, :pim]) {
      build_virtual_particle(:K_S0).by_closest_to(:K_S0)
      invariant_mass_of(:pip, :pim).within(0.487, 0.511)
      remove_used_particle_from_candidate_list
    }
   .kalman_kinematic_fit([:gamma, :gamma]) {
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
      chi2_cut 25
      neta ">=1"
    }
   .kinematic_fit([:prp, :K_S0, :eta]) {
      nominal
      constrain_four_momentum
      chi2_cut 200
    }

alg.with_decay_card(decay_card).apply(sel)
alg.execute_on([data_4600, incMC_4600, exMC_signal])