# ======================================================================
# Hybrid PBK model equations for nano/microplastic inhalation exposure
# Human PBK model: perfusion-limited and diffusion-limited submodels
#
# This script contains only the model equation structure used in the
# inhalation PBK framework. It does not include fitted parameter values,
# uncertainty analysis, QSS detection, particle-size distributions,
# calibration data, deposition coefficient tables, or unpublished inputs.
#
# WARNING:
# The example values are arbitrary placeholders for demonstrating code
# execution. They should not be interpreted as calibrated, fitted, or
# biologically realistic parameter values.
#
# Required package:
#   deSolve
#
# Author: Ira Wardani
# License: <choose license, e.g. MIT, GPL-3, CC-BY-4.0>
# ======================================================================

library(deSolve)

# ----------------------------------------------------------------------
# 1. Optional deposition function
# ----------------------------------------------------------------------
# In the full model, size-dependent deposition fractions are calculated
# using fitted deposition functions. For the public equation-only version,
# deposition fractions can be supplied directly by the user.
#
# fua  = fraction of inhaled particles deposited in upper airways
# ftra = fraction of inhaled particles deposited in tracheobronchial region
# fpul = fraction of inhaled particles deposited in pulmonary region
#
# The remaining inhaled particles are not initialized in the PBK system.

make_deposition_fractions <- function(fua, ftra, fpul) {
  fua  <- max(0, min(1, fua))
  ftra <- max(0, min(1, ftra))
  fpul <- max(0, min(1, fpul))
  
  s <- fua + ftra + fpul
  
  if (s > 1) {
    fua  <- fua / s
    ftra <- ftra / s
    fpul <- fpul / s
  }
  
  list(
    fua = fua,
    ftra = ftra,
    fpul = fpul
  )
}


# ----------------------------------------------------------------------
# 2. Initial conditions: inhalation exposure
# ----------------------------------------------------------------------
# Inhaled particles are initialized in the respiratory tract according to
# the deposition fractions:
#
#   Nua(0)  = dose_abs * fua
#   Ntra(0) = dose_abs * ftra
#   Npul(0) = dose_abs * fpul
#
# For the diffusion-limited model, the pulmonary deposited dose is placed
# in the pulmonary tissue subcompartment, NpulT.

state0_perf_inhalation <- function(dose_abs, fua, ftra, fpul) {
  dep <- make_deposition_fractions(fua, ftra, fpul)
  
  dep_dose <- dose_abs * (dep$fua + dep$ftra + dep$fpul)
  
  c(
    Na = 0,
    Nv = 0,
    
    Nua  = dose_abs * dep$fua,
    Ntra = dose_abs * dep$ftra,
    Npul = dose_abs * dep$fpul,
    
    Nlym = 0,
    
    Nstom = 0,
    Nsit_free = 0,
    Nsit_agg = 0,
    Nlit = 0,
    
    Ngi = 0,
    Nk = 0,
    Nli = 0,
    Ns = 0,
    Nh = 0,
    Nbr = 0,
    Nres = 0,
    
    Nmuc = 0,
    Nurine = 0,
    Nbile = 0,
    Nfeces = 0,
    
    Ntotal = dep_dose
  )
}


state0_diff_inhalation <- function(dose_abs, fua, ftra, fpul) {
  dep <- make_deposition_fractions(fua, ftra, fpul)
  
  dep_dose <- dose_abs * (dep$fua + dep$ftra + dep$fpul)
  
  c(
    NaB = 0,
    NaPC = 0,
    NvB = 0,
    NvPC = 0,
    
    Nlym = 0,
    
    Nua  = dose_abs * dep$fua,
    Ntra = dose_abs * dep$ftra,
    
    NpulT = dose_abs * dep$fpul,
    NpulB = 0,
    NpulPC = 0,
    
    Nstom = 0,
    Nsit_free = 0,
    Nsit_agg = 0,
    Nlit = 0,
    
    NgiT = 0,
    NgiB = 0,
    NgiPC = 0,
    
    NkT = 0,
    NkB = 0,
    NkPC = 0,
    
    NliT = 0,
    NliB = 0,
    NliPC = 0,
    
    NsT = 0,
    NsB = 0,
    NsPC = 0,
    
    NhT = 0,
    NhB = 0,
    NhPC = 0,
    
    NbrT = 0,
    NbrB = 0,
    NbrPC = 0,
    
    NresT = 0,
    NresB = 0,
    NresPC = 0,
    
    Nmuc = 0,
    Nurine = 0,
    Nbile = 0,
    Nfeces = 0,
    
    Ntotal = dep_dose
  )
}


# ----------------------------------------------------------------------
# 3. Perfusion-limited PBK model for inhalation exposure
# ----------------------------------------------------------------------

pbk_model_perf_inhalation <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # ------------------------------------------------------------------
    # Plasma binding / free fraction
    # ------------------------------------------------------------------
    Bmax <- Calb * N * MwPS / Mwalb * 1000
    fu   <- KD / (Bmax + KD)
    
    # ------------------------------------------------------------------
    # Concentrations
    # ------------------------------------------------------------------
    Ca    <- Na / Va
    Cv    <- Nv / Vv
    
    Cua   <- Nua / Vua
    Ctra  <- Ntra / Vtra
    Cvpul <- Npul / (Vpul * Ppul)
    Clym  <- Nlym / (Vlym * Plym)
    
    Cstom <- Nstom / Vstom
    Csit  <- (Nsit_free + Nsit_agg) / Vsit
    Clit  <- Nlit / Vlit
    
    Cvgi  <- Ngi  / (Vgi  * Pgi)
    Cvk   <- Nk   / (Vk   * Pk)
    Cvli  <- Nli  / (Vli  * Pli)
    Cvs   <- Ns   / (Vs   * Ps)
    Cvh   <- Nh   / (Vh   * Ph)
    Cvbr  <- Nbr  / (Vbr  * Pbr)
    Cvres <- Nres / (Vres * Pres)
    
    # ------------------------------------------------------------------
    # Cardiac output
    # ------------------------------------------------------------------
    QC <- Qk + Qli + Qh + Qbr + Qres + Qgi + Qs
    
    # ------------------------------------------------------------------
    # Blood and pulmonary exchange
    # ------------------------------------------------------------------
    dNa <- QC * (Cvpul - Ca) * fu
    
    dNv <- (Qk * Cvk +
              Qli * Cvli +
              Qh * Cvh +
              Qbr * Cvbr +
              Qres * Cvres) * fu +
      Nlym * Klym_ven -
      QC * Cv * fu
    
    dNpul <- QC * (Cv - Cvpul) * fu +
      Ntra * Ktra_pul -
      Npul * Kpul_tra -
      Npul * Kpul_ua -
      Npul * Kpul_lym
    
    # ------------------------------------------------------------------
    # Lymphatic transfer
    # ------------------------------------------------------------------
    dNlym <- Npul * Kpul_lym +
      Ngi  * Kgi_lym +
      Nk   * Kk_lym +
      Nli  * Kli_lym +
      Ns   * Ks_lym +
      Nh   * Kh_lym +
      Nbr  * Kbr_lym +
      Nres * Kres_lym -
      Nlym * Klym_ven
    
    # ------------------------------------------------------------------
    # Respiratory tract clearance
    # ------------------------------------------------------------------
    dNua <- Npul * Kpul_ua -
      Nua * Kua_tra -
      Nua * Kua_stom -
      Nua * Kmuc
    
    dNtra <- Nua * Kua_tra +
      Npul * Kpul_tra -
      Ntra * Ktra_stom -
      Ntra * Ktra_pul
    
    # ------------------------------------------------------------------
    # Swallowing and gastrointestinal transfer
    # ------------------------------------------------------------------
    dNstom <- Nua * Kua_stom +
      Ntra * Ktra_stom -
      Nstom * Kstom_sit
    
    dNsit_free <- Nstom * Kstom_sit +
      Nli * Kbile +
      Nsit_agg * Kdiss -
      Nsit_free * (Ksit_lit + Kagg + Ksit_gi)
    
    dNsit_agg <- Nsit_free * Kagg -
      Nsit_agg * (Kdiss + Ksit_agg_lit)
    
    dNlit <- Nsit_free * Ksit_lit +
      Nsit_agg * Ksit_agg_lit -
      Nlit * Kfec
    
    dNgi <- Qgi * (Ca - Cvgi) * fu +
      Nsit_free * Ksit_gi -
      Ngi * Kgi_lym
    
    # ------------------------------------------------------------------
    # Systemic organs
    # ------------------------------------------------------------------
    dNk <- Qk * (Ca - Cvk) * fu -
      Nk * Kurine -
      Nk * Kk_lym
    
    dNli <- Qli * (Ca - Cvli) * fu +
      Qs  * Cvs  * fu +
      Qgi * Cvgi * fu -
      Nli * Kbile -
      Nli * Kli_lym
    
    dNs <- Qs * (Ca - Cvs) * fu -
      Ns * Ks_lym
    
    dNh <- Qh * (Ca - Cvh) * fu -
      Nh * Kh_lym
    
    dNbr <- Qbr * (Ca - Cvbr) * fu -
      Nbr * Kbr_lym
    
    dNres <- Qres * (Ca - Cvres) * fu -
      Nres * Kres_lym
    
    # ------------------------------------------------------------------
    # Excretion sinks
    # ------------------------------------------------------------------
    dNmuc   <- Nua * Kmuc
    dNurine <- Nk  * Kurine
    dNbile  <- Nli * Kbile
    dNfeces <- Nlit * Kfec
    
    # ------------------------------------------------------------------
    # Mass-balance tracking
    # ------------------------------------------------------------------
    dNtotal <- dNa + dNv + dNua + dNtra + dNpul + dNlym +
      dNstom + dNsit_free + dNsit_agg + dNlit +
      dNgi + dNk + dNli + dNs + dNh + dNbr + dNres +
      dNmuc + dNurine + dNbile + dNfeces
    
    list(c(
      dNa,
      dNv,
      dNua,
      dNtra,
      dNpul,
      dNlym,
      dNstom,
      dNsit_free,
      dNsit_agg,
      dNlit,
      dNgi,
      dNk,
      dNli,
      dNs,
      dNh,
      dNbr,
      dNres,
      dNmuc,
      dNurine,
      dNbile,
      dNfeces,
      dNtotal
    ))
  })
}


# ----------------------------------------------------------------------
# 4. Diffusion-limited / hybrid PBK model for inhalation exposure
# ----------------------------------------------------------------------
# Compartments are separated into:
#   T  = tissue subcompartment
#   B  = blood or vascular subcompartment
#   PC = particle-retaining or phagocytic-cell subcompartment
#
# Uptake into particle-retaining compartments is represented using
# time-dependent saturable uptake functions:
#
#   Kup_i(t) = Kmax_i * t^n_i / (K50_i^n_i + t^n_i)
#
# All parameter values must be supplied through the parameter vector.

pbk_model_diff_inhalation <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # ------------------------------------------------------------------
    # Plasma binding / free fraction
    # ------------------------------------------------------------------
    Bmax <- Calb * N * MwPS / Mwalb * 1000
    fu   <- KD / (Bmax + KD)
    
    # ------------------------------------------------------------------
    # Cardiac output
    # ------------------------------------------------------------------
    QC <- Qk + Qli + Qh + Qbr + Qres + Qgi + Qs
    
    # ------------------------------------------------------------------
    # Blood concentrations
    # ------------------------------------------------------------------
    Ca <- NaB / Va
    Cv <- NvB / Vv
    
    # ------------------------------------------------------------------
    # Airway and pulmonary concentrations
    # ------------------------------------------------------------------
    Cua   <- Nua / Vua
    Ctra  <- Ntra / Vtra
    Cvpul <- NpulB / VpulB
    Ctpul <- NpulT / VpulT
    
    # ------------------------------------------------------------------
    # GI lumen concentrations
    # ------------------------------------------------------------------
    Cstom <- Nstom / Vstom
    Csit  <- (Nsit_free + Nsit_agg) / Vsit
    Clit  <- Nlit / Vlit
    
    # ------------------------------------------------------------------
    # Organ concentrations
    # ------------------------------------------------------------------
    Cvgi  <- NgiB  / VgiB
    Ctgi  <- NgiT  / VgiT
    
    Cvk   <- NkB   / VkB
    Ctk   <- NkT   / VkT
    
    Cvli  <- NliB  / VliB
    Ctli  <- NliT  / VliT
    
    Cvs   <- NsB   / VsB
    Cts   <- NsT   / VsT
    
    Cvh   <- NhB   / VhB
    Cth   <- NhT   / VhT
    
    Cvbr  <- NbrB  / VbrB
    Ctbr  <- NbrT  / VbrT
    
    Cvres <- NresB / VresB
    Ctres <- NresT / VresT
    
    # ------------------------------------------------------------------
    # Time-dependent uptake into particle-retaining compartments
    # ------------------------------------------------------------------
    Kup_blood <- (Kmax_blood * time^nblood) /
      (K50_blood^nblood + time^nblood)
    
    Kup_pul <- (Kmax_pul * time^npul) /
      (K50_pul^npul + time^npul)
    
    Kup_gi <- (Kmax_gi * time^ngi) /
      (K50_gi^ngi + time^ngi)
    
    Kup_k <- (Kmax_k * time^nk) /
      (K50_k^nk + time^nk)
    
    Kup_li <- (Kmax_li * time^nli) /
      (K50_li^nli + time^nli)
    
    Kup_s <- (Kmax_s * time^ns) /
      (K50_s^ns + time^ns)
    
    Kup_h <- (Kmax_h * time^nh) /
      (K50_h^nh + time^nh)
    
    Kup_br <- (Kmax_br * time^nbr) /
      (K50_br^nbr + time^nbr)
    
    Kup_res <- (Kmax_res * time^nres) /
      (K50_res^nres + time^nres)
    
    # ------------------------------------------------------------------
    # Blood
    # ------------------------------------------------------------------
    dNaB <- QC * (Cvpul - Ca) * fu -
      Kup_blood * NaB +
      Kres_blood * NaPC
    
    dNaPC <- Kup_blood * NaB -
      Kres_blood * NaPC
    
    dNvB <- (Qk * Cvk +
               Qli * Cvli +
               Qh * Cvh +
               Qbr * Cvbr +
               Qres * Cvres) * fu +
      Nlym * Klym_ven -
      QC * Cv * fu -
      Kup_blood * NvB +
      Kres_blood * NvPC
    
    dNvPC <- Kup_blood * NvB -
      Kres_blood * NvPC
    
    # ------------------------------------------------------------------
    # Lymph
    # ------------------------------------------------------------------
    dNlym <- NpulT * Kpul_lym +
      NgiT  * Kgi_lym +
      NkT   * Kk_lym +
      NliT  * Kli_lym +
      NsT   * Ks_lym +
      NhT   * Kh_lym +
      NbrT  * Kbr_lym +
      NresT * Kres_lym -
      Nlym  * Klym_ven
    
    # ------------------------------------------------------------------
    # Respiratory tract clearance
    # ------------------------------------------------------------------
    dNua <- NpulT * Kpul_ua -
      Nua * Kua_tra -
      Nua * Kua_stom -
      Nua * Kmuc
    
    dNtra <- Nua * Kua_tra +
      NpulPC * Kres_tra -
      Ntra * Ktra_stom -
      Ntra * Ktra_pul
    
    # ------------------------------------------------------------------
    # Pulmonary region
    # ------------------------------------------------------------------
    dNpulT <- Xpul * QC * (Cvpul * fu) -
      Ctpul / Ppul +
      Ntra * Ktra_pul -
      NpulPC * Kpul_tra -
      Kup_pul * NpulT +
      Kres_pul * NpulPC -
      NpulT * Kpul_ua -
      NpulT * Kpul_lym
    
    dNpulB <- QC * (Cv - Cvpul) * fu -
      Xpul * QC * (Cvpul * fu) +
      Ctpul / Ppul
    
    dNpulPC <- Kup_pul * NpulT -
      Kres_pul * NpulPC
    
    # ------------------------------------------------------------------
    # Swallowing and gastrointestinal transfer
    # ------------------------------------------------------------------
    dNstom <- Nua * Kua_stom +
      Ntra * Ktra_stom -
      Nstom * Kstom_sit
    
    dNsit_free <- Nstom * Kstom_sit +
      NliT * Kbile +
      Nsit_agg * Kdiss -
      Nsit_free * (Ksit_lit + Kagg + Ksit_gi)
    
    dNsit_agg <- Nsit_free * Kagg -
      Nsit_agg * (Kdiss + Ksit_agg_lit)
    
    dNlit <- Nsit_free * Ksit_lit +
      Nsit_agg * Ksit_agg_lit -
      Nlit * Kfec
    
    # ------------------------------------------------------------------
    # GI tissue
    # ------------------------------------------------------------------
    dNgiT <- Xgi * Qgi * (Cvgi * fu) -
      Ctgi / Pgi -
      Kup_gi * NgiT +
      Kres_gi * NgiPC -
      NgiT * Kgi_lym
    
    dNgiB <- Qgi * (Ca - Cvgi) * fu -
      Xgi * Qgi * (Cvgi * fu) +
      Ctgi / Pgi +
      Nsit_free * Ksit_gi
    
    dNgiPC <- Kup_gi * NgiT -
      Kres_gi * NgiPC
    
    # ------------------------------------------------------------------
    # Kidney
    # ------------------------------------------------------------------
    dNkT <- Xk * Qk * (Cvk * fu) -
      Ctk / Pk -
      Kup_k * NkT +
      Kres_k * NkPC -
      NkT * Kk_lym
    
    dNkB <- Qk * (Ca - Cvk) * fu -
      Xk * Qk * (Cvk * fu) +
      Ctk / Pk -
      NkB * Kurine
    
    dNkPC <- Kup_k * NkT -
      Kres_k * NkPC
    
    # ------------------------------------------------------------------
    # Liver
    # ------------------------------------------------------------------
    dNliT <- Xli * Qli * (Cvli * fu) -
      Ctli / Pli -
      Kup_li * NliT +
      Kres_li * NliPC -
      NliT * Kbile -
      NliT * Kli_lym
    
    dNliB <- Qli * (Ca - Cvli) * fu +
      Qs  * Cvs  * fu +
      Qgi * Cvgi * fu -
      Xli * Qli * (Cvli * fu) +
      Ctli / Pli
    
    dNliPC <- Kup_li * NliT -
      Kres_li * NliPC
    
    # ------------------------------------------------------------------
    # Spleen
    # ------------------------------------------------------------------
    dNsT <- Xs * Qs * (Cvs * fu) -
      Cts / Ps -
      Kup_s * NsT +
      Kres_s * NsPC -
      NsT * Ks_lym
    
    dNsB <- Qs * (Ca - Cvs) * fu -
      Xs * Qs * (Cvs * fu) +
      Cts / Ps
    
    dNsPC <- Kup_s * NsT -
      Kres_s * NsPC
    
    # ------------------------------------------------------------------
    # Heart
    # ------------------------------------------------------------------
    dNhT <- Xh * Qh * (Cvh * fu) -
      Cth / Ph -
      Kup_h * NhT +
      Kres_h * NhPC -
      NhT * Kh_lym
    
    dNhB <- Qh * (Ca - Cvh) * fu -
      Xh * Qh * (Cvh * fu) +
      Cth / Ph
    
    dNhPC <- Kup_h * NhT -
      Kres_h * NhPC
    
    # ------------------------------------------------------------------
    # Brain
    # ------------------------------------------------------------------
    dNbrT <- Xbr * Qbr * (Cvbr * fu) -
      Ctbr / Pbr -
      Kup_br * NbrT +
      Kres_br * NbrPC -
      NbrT * Kbr_lym
    
    dNbrB <- Qbr * (Ca - Cvbr) * fu -
      Xbr * Qbr * (Cvbr * fu) +
      Ctbr / Pbr
    
    dNbrPC <- Kup_br * NbrT -
      Kres_br * NbrPC
    
    # ------------------------------------------------------------------
    # Rest of body
    # ------------------------------------------------------------------
    dNresT <- Xres * Qres * (Cvres * fu) -
      Ctres / Pres -
      Kup_res * NresT +
      Kres_res * NresPC -
      NresT * Kres_lym
    
    dNresB <- Qres * (Ca - Cvres) * fu -
      Xres * Qres * (Cvres * fu) +
      Ctres / Pres
    
    dNresPC <- Kup_res * NresT -
      Kres_res * NresPC
    
    # ------------------------------------------------------------------
    # Excretion sinks
    # ------------------------------------------------------------------
    dNmuc   <- Nua  * Kmuc
    dNurine <- NkB  * Kurine
    dNbile  <- NliT * Kbile
    dNfeces <- Nlit * Kfec
    
    # ------------------------------------------------------------------
    # Mass-balance tracking
    # ------------------------------------------------------------------
    dNtotal <- dNaB + dNaPC + dNvB + dNvPC +
      dNlym + dNua + dNtra +
      dNpulT + dNpulB + dNpulPC +
      dNstom + dNsit_free + dNsit_agg + dNlit +
      dNgiT + dNgiB + dNgiPC +
      dNkT + dNkB + dNkPC +
      dNliT + dNliB + dNliPC +
      dNsT + dNsB + dNsPC +
      dNhT + dNhB + dNhPC +
      dNbrT + dNbrB + dNbrPC +
      dNresT + dNresB + dNresPC +
      dNmuc + dNurine + dNbile + dNfeces
    
    list(c(
      dNaB,
      dNaPC,
      dNvB,
      dNvPC,
      dNlym,
      dNua,
      dNtra,
      dNpulT,
      dNpulB,
      dNpulPC,
      dNstom,
      dNsit_free,
      dNsit_agg,
      dNlit,
      dNgiT,
      dNgiB,
      dNgiPC,
      dNkT,
      dNkB,
      dNkPC,
      dNliT,
      dNliB,
      dNliPC,
      dNsT,
      dNsB,
      dNsPC,
      dNhT,
      dNhB,
      dNhPC,
      dNbrT,
      dNbrB,
      dNbrPC,
      dNresT,
      dNresB,
      dNresPC,
      dNmuc,
      dNurine,
      dNbile,
      dNfeces,
      dNtotal
    ))
  })
}


# ----------------------------------------------------------------------
# 5. Minimal example of model execution
# ----------------------------------------------------------------------
# This section is intentionally schematic.
# Replace the dummy parameter values with user-defined or literature-based
# values before running the model.
#
# The parameter vector below contains placeholder values only.
# It includes both perfusion-model and diffusion-model parameters so that
# both example simulations can run.

run_example <- TRUE

if (run_example) {
  
  times <- seq(0, 24, by = 0.1)
  
  dose_abs <- 1e7
  
  # Example deposition fractions for inhalation.
  # In the full research workflow, these can be size-dependent.
  fua  <- 0.30
  ftra <- 0.20
  fpul <- 0.50
  
  # --------------------------------------------------------------------
  # Dummy parameter vector
  # --------------------------------------------------------------------
  # These values are placeholders for demonstration only.
  # Replace them with literature-based, calibrated, or user-defined values.
  
  parameters <- c(
    # Binding / free fraction
    Calb = 1,
    N = 1,
    MwPS = 1,
    Mwalb = 1,
    KD = 1,
    
    # Perfusion-model volumes
    Va = 1,
    Vv = 1,
    Vua = 1,
    Vtra = 1,
    Vpul = 1,
    Vlym = 1,
    Vstom = 1,
    Vsit = 1,
    Vlit = 1,
    Vgi = 1,
    Vk = 1,
    Vli = 1,
    Vs = 1,
    Vh = 1,
    Vbr = 1,
    Vres = 1,
    
    # Diffusion-model vascular/tissue volumes
    VpulT = 1,
    VpulB = 1,
    VgiT = 1,
    VgiB = 1,
    VkT = 1,
    VkB = 1,
    VliT = 1,
    VliB = 1,
    VsT = 1,
    VsB = 1,
    VhT = 1,
    VhB = 1,
    VbrT = 1,
    VbrB = 1,
    VresT = 1,
    VresB = 1,
    
    # Partition coefficients
    Ppul = 1,
    Plym = 1,
    Pgi = 1,
    Pk = 1,
    Pli = 1,
    Ps = 1,
    Ph = 1,
    Pbr = 1,
    Pres = 1,
    
    # Blood flows
    Qk = 1,
    Qli = 1,
    Qh = 1,
    Qbr = 1,
    Qres = 1,
    Qgi = 1,
    Qs = 1,
    
    # Tissue exchange fractions for diffusion model
    Xpul = 1,
    Xgi = 1,
    Xk = 1,
    Xli = 1,
    Xs = 1,
    Xh = 1,
    Xbr = 1,
    Xres = 1,
    
    # Airway, lymphatic, GI, and excretion rate constants
    Klym_ven = 0.01,
    Kpul_tra = 0.01,
    Kpul_ua = 0.01,
    Kpul_lym = 0.01,
    Kua_tra = 0.01,
    Kua_stom = 0.01,
    Kmuc = 0.01,
    Ktra_stom = 0.01,
    Ktra_pul = 0.01,
    
    Kstom_sit = 0.01,
    Kdiss = 0.01,
    Ksit_lit = 0.01,
    Kagg = 0.01,
    Ksit_gi = 0.01,
    Ksit_agg_lit = 0.01,
    Kfec = 0.01,
    
    Kgi_lym = 0.01,
    Kk_lym = 0.01,
    Kli_lym = 0.01,
    Ks_lym = 0.01,
    Kh_lym = 0.01,
    Kbr_lym = 0.01,
    Kres_lym = 0.01,
    
    Kurine = 0.01,
    Kbile = 0.01,
    
    # Reversibility / release from particle-retaining compartments
    Kres_blood = 0.01,
    Kres_pul = 0.01,
    Kres_tra = 0.01,
    Kres_gi = 0.01,
    Kres_k = 0.01,
    Kres_li = 0.01,
    Kres_s = 0.01,
    Kres_h = 0.01,
    Kres_br = 0.01,
    Kres_res = 0.01,
    
    # Time-dependent uptake parameters for diffusion model
    Kmax_blood = 0.01,
    K50_blood = 1,
    nblood = 1,
    
    Kmax_pul = 0.01,
    K50_pul = 1,
    npul = 1,
    
    Kmax_gi = 0.01,
    K50_gi = 1,
    ngi = 1,
    
    Kmax_k = 0.01,
    K50_k = 1,
    nk = 1,
    
    Kmax_li = 0.01,
    K50_li = 1,
    nli = 1,
    
    Kmax_s = 0.01,
    K50_s = 1,
    ns = 1,
    
    Kmax_h = 0.01,
    K50_h = 1,
    nh = 1,
    
    Kmax_br = 0.01,
    K50_br = 1,
    nbr = 1,
    
    Kmax_res = 0.01,
    K50_res = 1,
    nres = 1
  )
  
  # --------------------------------------------------------------------
  # 5.1 Inhalation: perfusion-limited model
  # --------------------------------------------------------------------
  
  out_inhalation_perf <- ode(
    y = state0_perf_inhalation(
      dose_abs = dose_abs,
      fua = fua,
      ftra = ftra,
      fpul = fpul
    ),
    times = times,
    func = pbk_model_perf_inhalation,
    parms = parameters,
    method = "lsoda"
  )
  
  out_inhalation_perf <- as.data.frame(out_inhalation_perf)
  
  cat("\nInhalation perfusion-limited model:\n")
  print(head(out_inhalation_perf))
  
  
  # --------------------------------------------------------------------
  # 5.2 Inhalation: diffusion-limited / hybrid model
  # --------------------------------------------------------------------
  
  out_inhalation_diff <- ode(
    y = state0_diff_inhalation(
      dose_abs = dose_abs,
      fua = fua,
      ftra = ftra,
      fpul = fpul
    ),
    times = times,
    func = pbk_model_diff_inhalation,
    parms = parameters,
    method = "lsoda"
  )
  
  out_inhalation_diff <- as.data.frame(out_inhalation_diff)
  
  cat("\nInhalation diffusion-limited / hybrid model:\n")
  print(head(out_inhalation_diff))
}
