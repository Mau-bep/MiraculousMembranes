# Literature notes: finite-range adhesion shell vs Deserno's zero-range wrapping diagram

Workstream: literature. Written incrementally (a usage limit may cut it off); sections are appended as they are done.
Convention for every entry: "READ" = full text opened and searched (pdftotext of the PDF that WebFetch saved, or the
HTML), "ABSTRACT/SNIPPET" = only an abstract page or a web-search snippet. Nothing below is cited from memory
unless it is explicitly tagged "UNVERIFIED (memory)" and then it must not be quoted as a reference.
Our claim under test: a finite-range adhesion shell (range delta) leaves the contact condition
w = (kappa/2) (curvature jump)^2 unchanged but adds a NEGATIVE effective contact-line tension
tau_eff = -T_p w sqrt(delta/Dc), so the phase lines of the (w~, sigma~) diagram move to larger w~ like sqrt(delta).

---

## 1. Deserno's papers (priority 1)

### 1.1 M. Deserno, "Elastic deformation of a fluid membrane upon colloid binding", Phys. Rev. E 69, 031903 (2004); arXiv:cond-mat/0303656
URL opened: https://arxiv.org/pdf/cond-mat/0303656 (and abs page https://arxiv.org/abs/cond-mat/0303656 appeared in search).
Status: READ the full text, but of **arXiv v1 (31 Mar 2003)**, not the published PRE version (the PRE text was not
accessible; page numbers/equation numbers below are those of v1 and may differ in PRE).

What it says that bears on us:
* **Zero-range / contact-energy model, stated as a modelling choice, not discussed quantitatively.** Sec. II: adhesion is "driven by a
  contact energy per unit area, w"; "Since the description will not aim at a microscopic understanding, continuum elasticity
  theory is taken as a basis" (Sec. II, before Eq. 1). There is NO parameter for the range of the adhesion potential
  anywhere in the paper and no discussion of how a finite range would modify the phase lines. The only statement about the
  validity limit is in footnote [38]: once the membrane deformation occurs on length scales comparable to the bilayer
  thickness ("typically a few nanometers") neither Helfrich elasticity nor the idealised surface description is appropriate and the
  bilayer structure would need to be modelled in more detail. (That footnote is made in the context of the z -> 2 neck, not
  of the adhesion potential.)
* **Contact condition.** Sec. III B, Eq. (14): a psi0dot = 1 - sqrt(w~), "the balance between adhesion energy and elastic membrane
  deformation results in a boundary condition on the contact curvature [17]; for curved substrates this becomes [18]", ref [18] =
  U. Seifert, R. Lipowsky, PRA 42, 4768 (1990). Deserno stresses that (14) holds only at the equilibrium state; in the energy
  E(z) calculation the contact curvature is fixed by asymptotic flatness at fixed z, and dE/dz = 0 reproduces (14) ("this would
  be one way to derive it"). This is exactly the identity dE/dz = (1 - psi0dot)^2 - w~ in our CONTEXT.md. He also notes that
  because the bending stiffness forbids a kink the "notion of a contact angle only remains meaningful in an asymptotic sense and is
  replaced by the concept of contact curvature" (Introduction, refs [17-20]).
* **Line tension is explicitly discussed and dismissed as a phenomenological device.** Sec. III (start): one could approximate
  E_free by a "phenomenological line energy", "However, neither the relation between the line tension constant and the membrane
  properties kappa and sigma would be known, nor is the implied dependency on the degree of wrapping -- namely ~ sin(alpha) --
  supported by more careful studies". Footnote [30]: Lipowsky & Doebereiner (EPL 43, 219 (1998), ref [14]) proposed a contact-line
  bending term ~ sqrt(z(2-z)) acting as an energy barrier; Deserno says for sigma = 0 such a term is "rigorously absent
  (due to its catenoid shape)" and for sigma > 0 the barrier "does not have the proposed form, since the bending energy
  is generally not localized at the rim". So Deserno's *elastic* line energy of the free membrane is not a line tension; it
  is an extended deformation. This is consistent with our finding that the elastic problem has no line tension at leading
  order and that tau_eff is purely a finite-range effect that the paper does not contain.
* **Phase lines.** Fig. 2 caption / Sec. III C: W (w~ = 4, continuous onset of partial wrapping, unchanged by E_free because E_free is of
  higher than linear order in z for small z), E (bold line, discontinuous partially wrapped <-> fully enveloped transition, i.e.
  equal total energy of the two minima of E~(z) in [0,2]), S1 and S2 (short dashed lines, spinodals belonging to E: S1 where the
  metastable partially wrapped state disappears on approaching full wrapping from below, i.e. binding/envelopment onset on increasing w~;
  S2 where the fully wrapped state loses its barrier on decreasing w~), plus the dotted line w~ = 4 + 2 sigma~ (E_total(z=2) = 0).
  Energy barrier in Fig. 3 at the E line; for the equal-energy criterion the minimum in z in [0,2] of Eq. (4) with E_free from the
  shooting solution. The numbers "sigma~ = 0.53 -> w~_E = 5.1" were NOT searched for in the text (Fig. 2 is graphical); we rely on CONTEXT.md.
* sigma~ = 0 limit: exact result (catenoid, E_free = 0): free for w~ < 4, fully wrapped for w~ > 4, no barrier. A finite range
  smears this into a crossover (our analysis); the paper has nothing on this.
* Also relevant: "sigma~ ~ 1 corresponds to a ~ lambda" with lambda = sqrt(kappa/sigma); for typical numbers kappa = 20 kT, sigma = 0.02 dyn/cm,
  lambda = 64 nm; biological particle radii 10 nm - few hundred nm. Hence for a = 30 nm and a "few nm" interaction range the dimensionless
  range s = delta/a is of order 0.1: not negligible against our sqrt(s) prediction.

### 1.2 M. Deserno and T. Bickel, "Wrapping of a spherical colloid by a fluid membrane", Europhys. Lett. 62, 767 (2003); arXiv:cond-mat/0212421
URL opened: https://arxiv.org/pdf/cond-mat/0212421 (abs page: https://arxiv.org/abs/cond-mat/0212421; journal page appeared in search:
https://epljournal.edpsciences.org/10.1209/epl/i2003-00438-4, not opened).
Status: READ (arXiv v1, 17 Dec 2002, 607 lines of extracted text, searched for range / line tension / thickness / Seifert / Lipowsky / simulation;
sections 1-3 read, rest skimmed through keyword search).
* Same model: "contact energy per unit area, w", no range parameter. They state the same dismissal of the line-tension ansatz
  ("neither the relation between the line tension constant and the membrane properties kappa and sigma would be known, nor is the implied
  dependency on the degree of wrapping correct") and propose to determine the exact membrane profile instead.
* sigma = 0: wc = 2 kappa/a^2 (i.e. w~ = 4), "no energy barrier", catenoid; for sigma > 0 the free -> partially wrapped transition is continuous at the same wc for any
  tension, partial -> full is discontinuous with a barrier "which surprisingly originates predominantly from tension".
* No comparison with simulations and no discussion of potential range in this paper either; it cites only Adv. Phys. 46, 13 (Seifert 1997) and
  Juelicher & Seifert 1994 for the shape equations.

Take-away for the study: both papers use a hard, zero-range contact energy and say nothing about finite-range corrections;
they do say that a naive line-tension term is not the right description of the *elastic* deformation. The finite-range tau_eff
of our analysis is a different object (an adhesion-layer, i.e. potential-range, effect) and so does not contradict them.

---

## 2. Finite-range potentials, contact potential, contact curvature (priority 2)

### 2.1 U. Seifert and R. Lipowsky, "Adhesion of vesicles", Phys. Rev. A 42, 4768 (1990)
URL opened: https://www.mpikg.mpg.de/rl/P/archive/067.pdf (author's archive copy at MPI-KG; OCR of a scan, so symbols are garbled but the
text is readable). Status: READ (intro, boundary-condition section, discussion, footnotes); the numerical shape catalogue not read in detail.
* **The contact potential is introduced as an explicit approximation to a finite-range potential.** Intro: "In order to have a bound state, the
  effective interaction potential must exhibit a minimum at a finite distance z_n. This potential range is typically of the order of a few nm", while the
  vesicle radius is 0.1-10 micrometres, so "we will ignore spatial variations on the scale of the potential range" and "replace the microscopic interaction
  potential for adhesion by an effective contact potential" (a footnote there cites Helfrich & Servuss, Nuovo Cimento D 3, 137 (1984) for a contact potential and
  E. Evans, Biophys. J. 48, 175 (1985) for a microscopic interaction potential; I did not open those). The paper therefore says nothing quantitative
  about corrections in (range / radius), except the one-line outlook "the influence of long-ranged adhesion potentials can be treated within the same
  theoretical framework".
* **Contact condition.** With bending rigidity "the contact potential no longer determines the contact angle (which is always equal to pi) but the contact
  curvature"; the boundary condition (their Eq. 2) comes from varying the contact point (a transversality condition, they cite Courant-Hilbert) and holds
  "irrespective of the chosen ensemble". Footnote 14: "For adhesion at a curved wall, the boundary condition (2) reads C1 = (2W/kappa)^(1/2) + C1^wall",
  i.e. W = (kappa/2)(C1 - C1^wall)^2: exactly our "w = (kappa/2) (curvature jump)^2" (their kappa is Helfrich's, energy kappa/2 (C1+C2)^2, same as Deserno's Eq. 1).
  Footnote 14 also says the same condition "applies if one considers a stretchable membrane with an elastic energy term (k/2)(DeltaA/A0)^2 or adds a shear energy"
  (so it is independent of the details of the other elastic terms, in the same spirit as our "any U(d,d')" first integral).
* Not found: any statement that a finite range changes the contact condition or produces a line tension.

### 2.2 R. Lipowsky and U. Seifert, "Adhesion of vesicles and membranes", Mol. Cryst. Liq. Cryst. 202, 17-25 (1991)
URL opened: https://www.mpikg.mpg.de/rl/P/archive/074.pdf. Status: converted to text and searched by keyword only (range / line tension / contact potential / effective).
Restates the contact-potential model, effective Young-Dupre angle W = S(1 - cos phi_eff), and focuses on thermal renormalisation of W (RG, Monte Carlo)
with potential widths of nm. No finite-range line tension found by keyword search. Low value for us.

### 2.3 M. Deserno, M. M. Mueller and J. Guven, "Contact lines for fluid surface adhesion", arXiv:cond-mat/0703019 (2007)
(journal reference believed to be Phys. Rev. E 76, 011605 (2007): from memory, UNVERIFIED; cite the arXiv id.)
URL opened: https://arxiv.org/pdf/cond-mat/0703019. Status: READ abstract and introduction, keyword-searched the rest.
* Derives the contact-line boundary conditions (Young-Dupre, contact curvature for Helfrich membranes) geometrically for general Hamiltonians, and states that
  "except for capillary phenomena, these boundary conditions are not the manifestation of a local force balance"; Hamiltonians with higher derivatives "notice changes
  in slope or even curvature". This is the correct general reference for why the contact condition is a transversality (variation of the contact point) condition and
  why it involves the curvature jump. It is also the cleanest citation for the statement that the contact condition is the sharp (zero-range) limit.
* Intro, on the validity of the contact energy: "In the majority of cases the spatial extension of the surface ... exceeds the range of interaction ... by a large
  amount" (van der Waals, hydrophobic, screened electrostatic forces "typically extend over several nanometers", vesicles/droplets microns to mm); "Under these
  conditions the interaction is well approximated by a contact energy". Again an assumption statement, no estimate of the correction.

### 2.4 R. Capovilla and J. Guven, "Geometry of lipid vesicle adhesion", arXiv:cond-mat/0203336 (2002; published in Phys. Rev. E, volume/page not checked)
URL opened: https://arxiv.org/pdf/cond-mat/0203336. Status: READ abstract + keyword search. Zero-range contact interaction ("interaction energy proportional to the area
of contact"), derives the boundary geometry (normal force balance, strong-bonding limit, curvature asymmetry). Nothing on finite range. Low value (background for the
contact condition for general substrates).

### 2.5 Leads found but NOT accessible (403 from the publisher), not read
* "Contact-line bending energy controls phospholipid vesicle adhesion", Proc. R. Soc. A 480 (2284), 20230545 (2024) (authors not retrieved),
  https://royalsocietypublishing.org/rspa/article/480/2284/20230545/101141/Contact-line-bending-energy-controls-phospholipid : HTTP 403. A web-search snippet (not
  verified) says it formulates "a competition between adhesive energy and the bending-related line tension at the edge of the contact region": potentially the most
  directly relevant recent paper on a bending-induced contact-line energy. TODO for Mau: open it by hand.
* "Adhesion of elastic membranes. Part I: a generalized Tabor parameter", Proc. R. Soc. A 482 (2330), 20250932,
  https://royalsocietypublishing.org/rspa/article/482/2330/20250932/479595/Adhesion-of-elastic-membranes-Part-I-a-generalized : HTTP 403. Snippet (unverified) says
  the contact curvature and the size of the rounded contact region follow from the ratio of bending rigidity and contact potential. A Tabor parameter is by construction
  a ratio (range of the interaction) / (elastic boundary-layer length), i.e. the analogue of our delta / ell_c = delta sqrt(w~)/a: if its content is as I expect, it is
  the natural place for a formal statement of when the zero-range limit holds. UNVERIFIED; do not cite until read.

---

## 3. Simulations compared with the zero-range wrapping theory (priority 3)

### 3.1 T. Ruiz-Herrero, E. Velasco, M. F. Hagan, "Mechanisms of budding of nanoscale particles through lipid bilayers", J. Phys. Chem. B 116, 9595 (2012); arXiv:1202.4691
URL opened: https://arxiv.org/pdf/1202.4691 (PDF text extracted and searched; also an ar5iv summary page, https://ar5iv.labs.arxiv.org/html/1202.4691, which I do NOT trust).
Status: READ the abstract, the elastic-theory set-up, the phase-diagram comparison paragraph and the discussion; not the dynamics section.
* Model: implicit-solvent Cooke-type lipids (WCA tails, head groups, sigma = 0.9 nm) and a smooth sphere attracting the heads with a **shifted Lennard-Jones
  potential cut off at r_ph + s with r_ph = 3.5 sigma**, i.e. an interaction range of a few nm (~3 nm by my reading of the cutoff; the exact definition of s was not
  read), particle radii R = 6-12 sigma (9-36 nm). So the dimensionless range delta/R is roughly 0.3-0.5: far from the zero-range limit.
* Elastic theory they compare with is **not Deserno's exact profile solution**: it is their own toroidal-rim ansatz "following Deserno & Gelbart [J. Phys. Chem. B 106, 5543 (2002)]"
  (cap of curvature 1/R + toroidal rim), with tensionless binodal eps* = 2 kappa/R^2 (our w~ = 4) quoted in their Fig. 13 discussion.
* Result: "the theory and simulations agree to within about 0.2 kBT" but "the theoretical binodal is below the computational results" (simulation needs MORE adhesion).
  The authors attribute this to the neglected configurational entropy of the lipids and to the theory's "infinitesimally thin membrane" assumption; they do not
  mention a line tension or the interaction range. Near the binodal they find **long-lived partially wrapped states** (metastable by free-energy calculations) that no
  elastic theory with a sharp contact predicts (neither Deserno's nor Zhang-Nguyen's, which they say have no metastable partial wrapping on an infinite tensionless membrane).
* Relevance: (i) sign agrees with our prediction (finite range -> binding at larger adhesion), (ii) the magnitude 0.2 kT is small and cannot discriminate sqrt(delta) from other
  explanations, (iii) the "stalled/long-lived partial wrapping near the transition" is qualitatively what a smoothed first-order transition (S1 and S2 closing up) would give,
  but they did not interpret it that way. This is circumstantial support, not a test.

### 3.2 A. H. Bahrami, M. Raatz, J. Agudo-Canalejo, R. Michel, E. M. Curtis, C. K. Hall, M. Gradzielski, R. Lipowsky, T. R. Weikl, "Wrapping of nanoparticles by membranes", Adv. Colloid Interface Sci. 208, 214-224 (2014)
URL opened: https://www.mpikg.mpg.de/rl/P/archive/397.pdf (author archive copy; full text extracted, Secs. 1-3 read in detail, rest keyword-searched). Status: READ (Sec. 3).
**This is the most relevant paper found for our question.** It is a review, and its central message is that the potential range rho of the particle-membrane
adhesion potential "crucially affects the wrapping process" (abstract and Sec. 1: for nanoparticles "the interaction range of the adhesion potential can be several percent of
the particle diameter", while "previous theoretical investigations have been largely focused on ... negligibly small" range).
* Sec. 3.1: for a contact potential (rho = 0) a tensionless membrane wraps discontinuously at u = 2 (u = rescaled adhesion energy; u = w~/2 in Deserno's notation since both
  scale w with kappa/R^2), bound part spherical, unbound part a catenoid with zero bending energy: this is Deserno's sigma~ = 0 limit (they cite Deserno & Bickel and Deserno 2004 as [29], [30]).
* Sec. 3.2 and Figs. 1-2 (computations of Raatz, Lipowsky, Weikl, ref [39] below; Morse potential V(d) = U(e^{-2d/rho} - 2e^{-d/rho}), min -U at d = 0, steeply repulsive for d < 0 = a soft core,
  NOT our symmetric cos^2 shell; membrane = Helfrich, tensionless, axisymmetric energy minimisation, example rho = 0.01 R): "In contrast to the discontinuous wrapping for rho = 0, the fraction of the wrapped particle area increases
  continuously with u for a finite potential range"; the process is "centered around the value u = 2"; "With decreasing potential range rho, the wrapping process becomes more abrupt and is finally discontinuous in the limit rho -> 0";
  "an increase in the potential range rho requires larger values of the rescaled adhesion energy u for full wrapping and subsequent membrane crossing".
* Mechanism, stated in words and identical to our item (ii): "The minimum energy E is negative for rho > 0 and decreases both with increasing u and increasing rho. The decrease of E with
  increasing potential range results from a favorable interplay of bending and adhesion energies in the contact region in which the membrane detaches from the particle ... the membrane
  already approaches the catenoidal shape of the unbound membrane segment with zero bending energy but still gains adhesion energy due to the finite potential range"; the contact region "becomes wider" with rho.
  I.e. a negative energy of order (adhesion gain in the detachment zone) at the contact circle, growing with rho: qualitatively exactly our negative tau_eff. They give no formula (no sqrt(rho) law and no line-tension language) in this review.
* Caveat found: for non-spherical particles (Sec. 4) the full-wrapping threshold u* DEcreases with rho (Fig. 4(b)), because the membrane can cut corners in regions of high particle curvature; for the sphere this does not apply. Also the transition is rendered
  continuous for the finite ranges they show (example rho = 0.01 R, a few percent of the particle size at most): i.e. in the tensionless case the smoothing of the sigma~ = 0 transition (our "no transition, only a crossover") is a known result of the finite-range literature, not a new claim.
* The same paper (and its PRL companion below) uses DMD/MC simulations of triangulated vesicles with a SQUARE-WELL vertex-particle interaction of cutoff d_c (see 3.3), i.e. the same kind of discretised finite-range adhesion as our code.

### 3.3 A. H. Bahrami, R. Lipowsky, T. R. Weikl, "Tubulation and Aggregation of Spherical Nanoparticles Adsorbed on Vesicles", Phys. Rev. Lett. 109, 188102 (2012)
URL opened: https://www.mpikg.mpg.de/rl/P/archive/373.pdf. Status: converted to text, READ the abstract and the model paragraph, keyword-searched the rest.
Triangulated-vesicle Monte Carlo/energy minimisation; "a short-ranged square-well interaction potential. A vertex is bound to a particle with energy U_A if the distance of the vertex from the particle surface is within a cutoff distance d_c."
Cites Deserno 2004 (ref 4) and Mueller-Deserno-Guven (ref 12). No statement about the shift of the sphere-wrapping transition with d_c found by keyword search (the paper is about tubes). Useful only as a precedent for a vertex-cutoff discretisation of the adhesion shell on a triangulated membrane.

### 3.4 Not opened (identified from search results or from the references of the above; do not cite without reading)
* (Moved to 3.5 below: Raatz et al. was opened after all.)
* "Physical mechanisms of nanoparticle-membrane interactions: A coarse-grained study", J. Chem. Phys. 164, 084906 (2026), https://pubs.aip.org/aip/jcp/article/164/8/084906/3381086 (search result only; bioRxiv 10.1101/2025.10.31.685850). Unknown content.
* "Adhesion and Aggregation of Spherical Nanoparticles on Lipid Membranes", https://pubmed.ncbi.nlm.nih.gov/33120231/ (search result only).
* Not found/not searched in the budget: Dasgupta-Auth-Gompper, Fosnaric et al., Yue-Zhang, Smith-Jasnow-Balazs, Saric-Cacciuto, Vacha-Frenkel: I have no verified information on whether they compare with Deserno's phase lines.

---

## 4. Discrete-membrane codes and adhesion to a sphere (priority 4)
Only a search was done (no full text read). Mem3DG: C. Zhu, C. T. Lee, P. Rangamani, "Mem3DG: Modeling Membrane Mechanochemical Dynamics in 3D using Discrete Differential Geometry", arXiv:2111.04460
(https://arxiv.org/pdf/2111.04460; Biophysical Reports 2022 per my memory, UNVERIFIED). ABSTRACT/SNIPPET only: discrete Helfrich-Canham-Evans Hamiltonian with forces derived via DDG; search snippet does not mention adhesion to a
bead or any required mesh resolution relative to a potential range. Nothing usable on "mesh resolution needed for a finite-range adhesion shell". Ramakrishnan et al. not searched.
Only the generic point can be made from the sources above: a vertex-based square-well cutoff (Bahrami et al. PRL 2012) is the standard discretisation on triangulated membranes, and the continuum analyses that exist (Bahrami review Sec. 3.2, rho = 0.01 R in the example) take the range much SMALLER than the particle but then need a mesh well below rho; no paper found states a resolution criterion.

## 3.5 (added late) M. Raatz, R. Lipowsky, T. R. Weikl, "Cooperative wrapping of nanoparticles by membrane tubes", Soft Matter 10, 3570-3577 (2014), doi 10.1039/c3sm52498a; arXiv:1401.1668
URL opened: https://arxiv.org/pdf/1401.1668 (PDF text extracted and searched; Secs. on the single particle read, tube parts only skimmed). Status: READ the single-particle part.
This is the source (ref [39]) of the finite-range calculations quoted in the review of 3.2.
* Method: Helfrich energy, tensionless membrane, rotationally symmetric profiles r(z), z(r) discretised with "up to 1000 discretization points", derivatives by finite differences, constrained
  minimisation in Mathematica; Morse adhesion potential of range rho (examples rho/R = 0.01, 0.1, ...), u = U R^2/kappa.
* Single sphere: "The minimum energy E is negative for rho > 0 and decreases with u" and "also decreases with increasing potential range rho because the interplay of bending and adhesion that leads to the
  local minimum in the energy profiles e(z) [in the contact region] is more pronounced for larger values of rho"; "In the limit rho -> 0, the total energy E tends towards ... E = 0 for u < 2 and E = 4 pi kappa (2-u) for u > 2" (Deserno's sigma~ = 0 result). Wrapping "continuously increases with u", "centered around u = 2", "more abrupt" for smaller rho.
* **Directly the object we want**: at u = 2 (our w~ = 4) the energy density is zero both on the bound sphere and on the unbound catenoid, so the whole energy at u = 2 is the contact-region energy, and it is negative and "the minima in the energy
  profiles become broader" with rho (Sec. on energy densities, Fig. 6). That is the excess energy of the contact line in our language: E(u = 2; rho) < 0, with the width of the contact region growing with rho.
  The paper gives numbers only graphically (Fig. 3(c)), no scaling law in rho, no "line tension" language. TODO for the team: digitise E(u = 2) vs rho/R from Fig. 3(c) and compare with tau_eff = -T_p w sqrt(delta/Dc) (note the Morse potential has a soft inner repulsion, so T_p differs from our symmetric cos^p shell, and delta must be mapped through the Morse decay length).

## 2.6 Lead worth reading (not opened): U. Seifert, "Self-consistent theory of bound vesicles", Phys. Rev. Lett. 74, 5060 (1995)
Cited by Deserno 2004 as ref [20] for the contact-curvature concept (verified from his reference list). A web search returned only the citation (no abstract available). From memory (UNVERIFIED) it derives the contact curvature
from a microscopic potential self-consistently, which would be the natural place for a statement like "the contact condition is renormalised by O(range/radius)". Not read. Also unread: R. Lipowsky, U. Seifert, Langmuir 7, 1867 (1991) (cited by Deserno as [19]) and the Handbook chapter by Seifert & Lipowsky (1995), both containing the finite-range discussion in the Seifert-Lipowsky tradition.

---

## 5. Synthesis

**Does the literature support "finite range = negative line tension that shifts the lines like sqrt(range)"?** Qualitatively yes, quantitatively not yet. What is directly supported by sources I opened: (i) the zero-range contact condition w = (kappa/2)(C1 - C1^wall)^2 and its status as the
limit of a finite-range potential much shorter than every elastic scale (Seifert & Lipowsky 1990, incl. their footnote on a curved wall; Deserno-Mueller-Guven 2007; Deserno 2004 Eq. 14 and the dE/dz = 0 derivation); (ii) for a tensionless membrane and a finite-range (Morse) potential,
full Helfrich minimisation (Raatz, Lipowsky, Weikl 2014; reviewed by Bahrami et al. 2014) shows exactly the three features of our analysis: the minimum energy is negative for rho > 0 because in the detachment region the membrane is already catenoidal (zero bending) but still gains adhesion
(mechanism (ii) of our CONTEXT bookkeeping, the negative contact-line energy; the energy at u = 2 equals this contact-region energy and grows in magnitude with the width of the contact region), the sigma~ = 0 discontinuous transition becomes a continuous crossover centred at w~ = 4 (u = 2), more abrupt the smaller the range (our item 5 of line_tension_analysis.md,
and the reason the "S1 = 5.4 vs 4 at sigma~ = 0" is a definition artefact), and a larger range requires larger adhesion for full wrapping (our shift of E, S1, S2 to larger w~); (iii) the only direct MD comparison with an elastic theory that I could read (Ruiz-Herrero, Velasco, Hagan 2012, range ~ 3 nm on 9-36 nm particles)
finds the simulated wrapping binodal above the theoretical one by about 0.2 kT, i.e. the same sign, though the authors attribute it to lipid entropy and membrane thickness and see long-lived partially wrapped states near the transition. What I did NOT find in any opened source: the words "line tension" tied to the potential range, a sqrt(range) law,
the statement that the contact condition is unchanged for a finite-range shell at leading order, or a quantitative comparison with Deserno's sigma~ > 0 phase lines (E, S1, S2). Deserno 2004 never discusses the range and, in Sec. III and footnote [30], argues that a phenomenological line tension is not the right description of the *elastic* free-membrane energy (which for sigma = 0 is exactly zero for a catenoid, and for sigma > 0 is not
localised at the rim): this is consistent with, but different from, our tau_eff, which is a finite-range (adhesion-layer) effect. So the sqrt(range) scaling, the unchanged contact condition and the sign/size of tau_eff are, as far as this search can tell, results of this study that need their own derivation (our first integral and boundary-layer scaling) and a numerical check.
**What to cite for a formal defence:** for the zero-range limit and the contact condition: Seifert & Lipowsky, PRA 42, 4768 (1990) [read], Deserno, PRE 69, 031903 (2004) [arXiv v1 read], Deserno-Mueller-Guven, arXiv:cond-mat/0703019 [intro read]; for the reference diagram: Deserno & Bickel, EPL 62, 767 (2003); for external evidence that finite range makes the energy negative, smooths the sigma~ = 0 transition and moves full wrapping to larger u:
Raatz-Lipowsky-Weikl, Soft Matter 10, 3570 (2014) and Bahrami et al., Adv. Colloid Interface Sci. 208, 214 (2014), Sec. 3.2; for an MD precedent: Ruiz-Herrero-Velasco-Hagan, J. Phys. Chem. B 116, 9595 (2012); for the discretisation precedent (vertex cutoff on a triangulated membrane): Bahrami-Lipowsky-Weikl, PRL 109, 188102 (2012). The cleanest quantitative external test available is to digitise Fig. 3(c) of Raatz et al. (E at u = 2 versus rho/R) and compare with tau_eff.
**Could not access / not done:** the published PRE version of Deserno 2004 (only arXiv v1 read, equation numbers may differ); the Proc. R. Soc. A papers on contact-line bending energy and the generalized Tabor parameter (HTTP 403; potentially the most relevant recent literature, to be opened by hand); Seifert PRL 74, 5060 (1995), Lipowsky-Seifert Langmuir 1991, the Lipowsky/Seifert Handbook chapter; wetting/line-tension literature on finite-range interfacial potentials (not searched,
the natural analogy for "range renormalises into a line tension"); simulation papers by Dasgupta-Auth-Gompper, Fosnaric et al., Yue-Zhang, Smith-Jasnow-Balazs, Saric-Cacciuto, Vacha-Frenkel (not searched or only snippets); Mem3DG / Ramakrishnan et al. (abstract-level only, nothing on resolution requirements); numerical values from Fig. 2 of Deserno (graphical). I exceeded the nominal ~40-call budget by roughly a dozen calls to open the Raatz/Bahrami papers because they were the best hit.
