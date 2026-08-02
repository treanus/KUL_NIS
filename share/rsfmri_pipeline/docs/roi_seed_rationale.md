# ROI and Seed Selection Rationale
## Hierarchical rsfMRI Analysis Pipeline

This document describes the anatomical seeds and resting-state networks selected
for each condition profile, with the methodological justification and supporting
literature for each choice.

---

## 1. Foundational Resting-State Networks

Before condition-specific profiles, the pipeline defines a set of canonical
resting-state networks (RSNs) whose seed definitions are shared across profiles.
These were established in the seminal work of Biswal et al. (1995), who first
demonstrated that spontaneous BOLD fluctuations in the motor cortex exhibit
coherent low-frequency structure at rest, and consolidated by large-scale
decompositions such as Smith et al. (2009) and the Yeo et al. (2011)
7- and 17-network parcellations.

### 1.1 Sensorimotor Network (SMN)

**Seed:** Bilateral precentral and postcentral gyri (Harvard-Oxford Cortical atlas).

The SMN reflects the somatotopic organisation of primary motor (M1) and
somatosensory (S1) cortices — the cortical homunculus described by Penfield &
Rasmussen (1950). At rest, motor cortex BOLD fluctuations are highly coherent
bilaterally, forming one of the most robust and replicable RSNs (Biswal et al.,
1995). Within-mask ICA is expected to resolve 5 subcomponents corresponding to
the major somatotopic subdivisions: bilateral foot, hand, face/lips, trunk, and
supplementary motor area (SMA). This decomposition has been consistently observed
in both task-based and resting-state studies (Yeo et al., 2011; Deco et al., 2011).

**References:**
- Biswal B, Yetkin FZ, Haughton VM, Hyde JS (1995). Functional connectivity in
  the motor cortex of resting human brain using echo-planar MRI. *Magnetic
  Resonance in Medicine*, 34(4), 537–541.
- Penfield W, Rasmussen T (1950). *The Cerebral Cortex of Man*. Macmillan.
- Yeo BTT et al. (2011). The organization of the human cerebral cortex estimated
  by intrinsic functional connectivity. *Journal of Neurophysiology*, 106(3),
  1125–1165.

---

### 1.2 Default Mode Network (DMN)

**Seed:** Posterior cingulate cortex / precuneus + medial prefrontal cortex
(Harvard-Oxford Cortical atlas).

The DMN was defined by Raichle et al. (2001) as a set of regions showing
consistent deactivation during externally directed cognitive tasks, and later
characterised as a coherent resting-state network by Greicius et al. (2003).
The PCC/precuneus is used as the primary seed because it is the most reliably
connected hub of the DMN across individuals and scanners (Andrews-Hanna et al.,
2010; Buckner et al., 2008). Within-mask ICA is expected to resolve 5
subcomponents corresponding to the established DMN subsystems: anterior (mPFC),
core posterior (PCC/precuneus), lateral temporal, medial temporal lobe (MTL/
hippocampal), and angular gyrus subsystems (Andrews-Hanna et al., 2010).

**References:**
- Raichle ME et al. (2001). A default mode of brain function. *PNAS*, 98(2),
  676–682.
- Greicius MD, Krasnow B, Reiss AL, Menon V (2003). Functional connectivity in
  the resting brain: a network analysis of the default mode hypothesis. *PNAS*,
  100(1), 253–258.
- Buckner RL, Andrews-Hanna JR, Schacter DL (2008). The brain's default network:
  anatomy, function, and relevance to disease. *Annals of the New York Academy
  of Sciences*, 1124, 1–38.
- Andrews-Hanna JR, Reidler JS, Sepulcre J, Poulin R, Buckner RL (2010).
  Functional-anatomic fractionation of the brain's default network. *Neuron*,
  65(4), 550–562.

---

### 1.3 Dorsal Attention Network (DAN)

**Seed:** Bilateral superior parietal lobule / IPS (Harvard-Oxford Cortical atlas).

The DAN — comprising the frontal eye fields (FEF) and intraparietal sulcus (IPS)
— supports top-down, goal-directed attentional control (Corbetta & Shulman, 2002).
It is typically anticorrelated with the DMN at rest (Fox et al., 2005). The IPS
is used as seed because it is the more anatomically stable anchor across
individuals; FEF is typically recoverable as a subcomponent of within-mask ICA.

**References:**
- Corbetta M, Shulman GL (2002). Control of goal-directed and stimulus-driven
  attention in the brain. *Nature Reviews Neuroscience*, 3(3), 201–215.
- Fox MD et al. (2005). The human brain is intrinsically organized into dynamic,
  anticorrelated functional networks. *PNAS*, 102(27), 9673–9678.

---

### 1.4 Language Network

**Seed:** Left IFG pars triangularis and pars opercularis / Broca's area
(BIP Language Atlas; Harvard-Oxford Cortical atlas as fallback).

Language processing involves a left-lateralised network spanning inferior frontal
(Broca's area, BA44/45) and posterior temporal (Wernicke's area, BA22/STG/STS)
cortices, connected via the arcuate fasciculus (Hickok & Poeppel, 2007). The
dual-stream model (Price, 2010) distinguishes a dorsal stream (articulatory/
syntactic, IFG–premotor) from a ventral stream (semantic, temporal pole–IFG).
The BIP Language Atlas (Duffau et al.) provides a more surgical-grade
parcellation of language territories than generic cortical atlases, making it
preferable for pre-surgical language mapping. Broca's area is used as the primary
seed given its role as the key production hub; Wernicke's area typically emerges
as a subcomponent in within-mask ICA.

**References:**
- Hickok G, Poeppel D (2007). The cortical organization of speech processing.
  *Nature Reviews Neuroscience*, 8(5), 393–402.
- Price CJ (2010). The anatomy of language: a review of 100 fMRI studies
  published in 2009. *Annals of the New York Academy of Sciences*, 1156, 54–97.
- Duffau H et al. (2014). Intraoperative mapping of the subcortical language
  pathways using direct stimulations. *Brain*, 137(3), 828–841.

---

### 1.5 Salience Network

**Seed:** Bilateral anterior insula + dorsal anterior cingulate cortex
(Harvard-Oxford Cortical atlas).

The Salience Network (SN), anchored in the anterior insula and dACC, mediates
detection of and response to behaviourally relevant stimuli, and acts as a switch
between the DMN and task-positive networks (Menon & Uddin, 2010; Seeley et al.,
2007). The anterior insula seed captures the visceromotor and interoceptive
functions central to this network.

**References:**
- Seeley WW et al. (2007). Dissociable intrinsic connectivity networks for
  salience processing and executive control. *Journal of Neuroscience*, 27(9),
  2349–2356.
- Menon V, Uddin LQ (2010). Saliency, switching, attention and control: a
  network model of insula function. *Brain Structure and Function*, 214(5–6),
  655–667.

---

### 1.6 Visual Network

**Seed:** Bilateral intracalcarine cortex / lingual gyrus (Harvard-Oxford
Cortical atlas).

Primary visual cortex (V1/V2) in and around the calcarine sulcus provides the
most stable anchor for the visual RSN. Within-mask ICA is expected to resolve
dorsal stream (superior occipital / parietal), ventral stream (fusiform /
inferior temporal), and MT+ motion components.

---

### 1.7 Mesolimbic / Reward Network

**Seed:** Nucleus accumbens (HO Subcortical + FSL Striatum limbic subdivision),
amygdala, hippocampus (HO Subcortical).

The mesolimbic system — comprising the ventral tegmental area (VTA), NAcc,
amygdala, hippocampus, and OFC — underpins reward processing, motivational
salience, and emotional memory (Haber & Knutson, 2010). The NAcc is the
canonical functional seed, with the FSL Striatum limbic subdivision providing
connectivity-informed boundaries that more accurately capture the ventral
striatum than anatomical borders alone.

**References:**
- Haber SN, Knutson B (2010). The reward circuit: linking primate anatomy and
  human imaging. *Neuropsychopharmacology*, 35(1), 4–26.

---

## 2. Condition-Specific Profiles

---

### 2.1 Anhedonia / Reward Deficit

**Seeds:** NAcc (FSL Striatum limbic + HO Subcortical), amygdala, hippocampus,
OFC (HO Cortical).
**Networks:** Mesolimbic, DMN, Salience.

Anhedonia — the diminished capacity to experience pleasure — is a transdiagnostic
symptom central to MDD, schizophrenia, and substance use disorders. It is
mechanistically linked to disrupted signalling in the mesolimbic dopaminergic
circuit, particularly the VTA–NAcc pathway (Treadway & Zald, 2011). Resting-state
studies consistently report reduced NAcc connectivity with the prefrontal cortex
and altered amygdala–mPFC coupling in anhedonic patients (Pizzagalli, 2014).
The OFC seed is included because reward valuation and effort-cost computation
depend critically on OFC–striatum communication (Haber & Knutson, 2010).

**References:**
- Treadway MT, Zald DH (2011). Reconsidering anhedonia in depression: lessons
  from translational neuroscience. *Neuroscience & Biobehavioral Reviews*, 35(3),
  537–555.
- Pizzagalli DA (2014). Depression, stress, and anhedonia: toward a synthesis and
  integrated model. *Annual Review of Clinical Psychology*, 10, 393–423.

---

### 2.2 Major Depressive Disorder (MDD)

**Seeds:** Subgenual ACC / sgACC (HO Cortical anterior cingulate), amygdala,
hippocampus, NAcc (HO Subcortical + FSL Striatum limbic), DLPFC (HO Cortical
middle frontal), OFC (HO Cortical frontal orbital).
**Networks:** DMN, Mesolimbic, Salience, Frontoparietal.

The sgACC (Brodmann area 25) is the defining seed for MDD. It is hyperactive in
depression (Drevets et al., 1997), shows normalisation with successful treatment
across modalities (antidepressants, CBT, ECT), and is the direct target of
deep brain stimulation for treatment-resistant depression (Mayberg et al., 2005)
— one of the landmark findings in modern psychiatry. Pathological sgACC
hyperconnectivity with the DMN and disrupted connectivity with DLPFC are among
the most replicated functional findings in MDD (Greicius et al., 2007;
Hamilton et al., 2011). The amygdala seed captures the sustained emotional
reactivity characteristic of depression; hippocampal atrophy (driven by
stress/glucocorticoid exposure) motivates its inclusion as a connectivity seed
(Campbell & Macqueen, 2004). DLPFC is included as the primary target of
repetitive TMS in MDD, whose therapeutic effect is thought to be mediated via
its anticorrelated connectivity with the sgACC (Fox et al., 2012).

**References:**
- Mayberg HS et al. (2005). Deep brain stimulation for treatment-resistant
  depression. *Neuron*, 45(5), 651–660.
- Drevets WC et al. (1997). Subgenual prefrontal cortex abnormalities in mood
  disorders. *Nature*, 386(6627), 824–827.
- Greicius MD et al. (2007). Resting-state functional connectivity in major
  depression: abnormally increased contributions from subgenual cingulate cortex
  and thalamus. *Biological Psychiatry*, 62(5), 429–437.
- Hamilton JP et al. (2011). Functional neuroimaging of major depressive disorder:
  a meta-analysis and new integration of baseline activation and neural response
  data. *American Journal of Psychiatry*, 168(2), 114–122.
- Fox MD et al. (2012). Clinical applications of resting state fMRI: using
  inter-region correlations to target neuromodulatory treatments. *Dialogues in
  Clinical Neuroscience*, 14(4), 385–396.
- Campbell S, Macqueen G (2004). The role of the hippocampus in the pathophysiology
  of major depression. *Journal of Psychiatry & Neuroscience*, 29(6), 417–426.

---

### 2.3 Bipolar Disorder

**Seeds:** Amygdala (HO Subcortical), VLPFC / IFG (HO Cortical), OFC (HO
Cortical), ACC — dorsal and subgenual (HO Cortical), thalamus pre-frontal
subdivision (FSL Thalamus), hippocampus (HO Subcortical).
**Networks:** DMN, Salience, Limbic, Frontoparietal.

Bipolar disorder involves dysregulation of frontolimbic circuits mediating
emotion regulation, with the amygdala and ventrolateral PFC as primary nodes
(Phillips & Swartz, 2014). The amygdala is hyperactive and shows reduced
top-down inhibitory connectivity with the VLPFC in both manic and depressive
phases (Rich et al., 2006). Thalamic dysconnectivity — particularly involving
prefrontal thalamic nuclei — is a distinguishing feature of bipolar relative to
MDD and schizophrenia (Anticevic et al., 2014), making the FSL Thalamus atlas
pre-frontal subdivision an important differentiating seed. Hippocampal volume
loss and connectivity disruption parallel those seen in MDD but may reflect
distinct pathophysiology involving lithium-sensitive signalling cascades
(Strakowski et al., 2012). The inclusion of dorsal as well as subgenual ACC
seeds distinguishes bipolar from MDD, where sgACC involvement is more selective.

**References:**
- Phillips ML, Swartz HA (2014). A critical appraisal of neuroimaging studies of
  bipolar disorder: toward a new conceptualization of underlying neural circuitry
  and roadmap for future research. *American Journal of Psychiatry*, 171(8),
  829–843.
- Rich BA et al. (2006). Limbic hyperactivation during processing of neutral
  facial expressions in children with bipolar disorder. *PNAS*, 103(23),
  8900–8905.
- Anticevic A et al. (2014). Characterizing thalamo-cortical disturbances in
  schizophrenia and bipolar illness. *Cerebral Cortex*, 24(12), 3116–3130.
- Strakowski SM, Adler CM, Almeida J, Altshuler LL, Blumberg HP, Chang KD et al.
  (2012). The functional neuroanatomy of bipolar disorder: a consensus model.
  *Bipolar Disorders*, 14(4), 313–325.

---

### 2.4 ADHD

**Seeds:** Caudate executive subdivision (FSL Striatum), putamen rostral-motor
(FSL Striatum), NAcc/limbic striatum (FSL Striatum + HO Subcortical), DLPFC
(HO Cortical middle frontal), IFG (HO Cortical, for response inhibition), ACC
(HO Cortical), cerebellum (SUIT atlas — to be added).
**Networks:** Frontostriatal, DAN, DMN, Salience.

ADHD is characterised by dysfunction in two partially dissociable circuits:
a frontostriatal circuit subserving executive control (caudate–DLPFC), and a
mesolimbic circuit subserving reward and motivation (NAcc–OFC) (Castellanos &
Proal, 2012). The dorsal caudate and DLPFC are the canonical nodes of the
executive loop, consistently showing reduced connectivity in ADHD (Dickstein
et al., 2006). Response inhibition deficits implicate the IFG (right-lateralised
stop-signal network; Aron & Poldrack, 2006). A second major pathophysiological
model proposes that spontaneous DMN activity abnormally intrudes into task
processing (the "default mode interference" hypothesis), explaining
inattentiveness as a failure to suppress DMN at task onset (Sonuga-Barke &
Castellanos, 2007). Large-scale meta-analyses confirm reduced DAN and increased
DMN connectivity at rest (Cortese et al., 2012). Cerebellar-cortical connectivity
is disrupted in ADHD, consistent with timing and motor coordination deficits
(Stoodley, 2014), motivating addition of a SUIT-based cerebellar seed when that
atlas is integrated.

**References:**
- Castellanos FX, Proal E (2012). Large-scale brain systems in ADHD: beyond the
  prefrontal-striatal model. *Trends in Cognitive Sciences*, 16(1), 17–26.
- Dickstein SG, Bannon K, Castellanos FX, Milham MP (2006). The neural correlates
  of attention deficit hyperactivity disorder: an ALE meta-analysis. *Journal of
  Child Psychology and Psychiatry*, 47(10), 1051–1062.
- Sonuga-Barke EJS, Castellanos FX (2007). Spontaneous attentional fluctuations
  in impaired states and pathological conditions: a neurobiological hypothesis.
  *Neuroscience & Biobehavioral Reviews*, 31(7), 977–986.
- Cortese S et al. (2012). Toward systems neuroscience of ADHD: a meta-analysis
  of 55 fMRI studies. *American Journal of Psychiatry*, 169(10), 1038–1055.
- Aron AR, Poldrack RA (2006). Cortical and subcortical contributions to stop
  signal response inhibition: role of the subthalamic nucleus. *Journal of
  Neuroscience*, 26(9), 2424–2433.
- Stoodley CJ (2014). Distinct regions of the cerebellum show gray matter
  decreases in autism, ADHD, and developmental dyslexia. *Frontiers in Systems
  Neuroscience*, 8, 92.

---

### 2.5 OCD

**Seeds:** Caudate executive subdivision (FSL Striatum), OFC (HO Cortical
frontal orbital), ACC dorsal (HO Cortical), DLPFC (HO Cortical middle frontal),
thalamus pre-frontal subdivision (FSL Thalamus).
**Networks:** Frontostriatal, Salience, Frontoparietal.

OCD is the prototypical disorder of the cortico-striato-thalamo-cortical (CSTC)
loop (Saxena & Rauch, 2000). The direct pathway (striatum → GPi/SNr → thalamus →
cortex) is hyperactive, producing intrusive thoughts and compulsive behaviours;
the indirect pathway is relatively underactive, failing to suppress inappropriate
actions. The orbitofronto-caudate circuit is the most consistently implicated
(Rotge et al., 2008; Milad & Rauch, 2012), with OFC hyperactivity correlating
with symptom severity and normalising with successful treatment (CBT or SRIs).
The thalamus serves as the output relay of the CSTC loop and shows structural
and functional abnormalities in OCD (Anticevic et al., 2014). DLPFC is included
as the seat of cognitive control deficits and a target for TMS augmentation.

**References:**
- Saxena S, Rauch SL (2000). Functional neuroimaging and the neuroanatomy of
  obsessive-compulsive disorder. *Psychiatric Clinics of North America*, 23(3),
  563–586.
- Rotge JY et al. (2008). Provocation of obsessive-compulsive symptoms: a
  quantitative voxel-based meta-analysis of functional neuroimaging studies.
  *Journal of Psychiatry & Neuroscience*, 33(5), 405–412.
- Milad MR, Rauch SL (2012). Obsessive-compulsive disorder: beyond segregated
  cortico-striatal pathways. *Trends in Cognitive Sciences*, 16(1), 43–51.

---

### 2.6 ASD (Autism Spectrum Disorder)

**Seeds:** STS / superior temporal sulcus (HO Cortical), TPJ / temporoparietal
junction (HO Cortical angular gyrus + supramarginal gyrus), anterior insula
(HO Cortical), PCC (HO Cortical, DMN anchor).
**Networks:** DMN, Salience, Social brain.

The "underconnectivity theory" of ASD proposes reduced long-range functional
connectivity, particularly within networks subserving social cognition (Just et
al., 2004). The STS is a critical hub for biological motion processing, theory
of mind, and voice recognition — all impaired in ASD — and consistently shows
reduced activation and connectivity (Pelphrey et al., 2005). The TPJ, a node
shared by the mentalising and attention networks, shows atypical functional
lateralisation and connectivity in ASD (Lombardo et al., 2011). DMN connectivity
is disrupted, with reduced PCC–mPFC coherence correlating with social impairment
severity (Kennedy & Courchesne, 2008). Anterior insula is included given its
role in interoception and its altered connectivity within the salience network
in ASD (Uddin & Menon, 2009).

**References:**
- Just MA, Cherkassky VL, Keller TA, Minshew NJ (2004). Cortical activation and
  synchronization during sentence comprehension in high-functioning autism: evidence
  of underconnectivity. *Brain*, 127(8), 1811–1821.
- Kennedy DP, Courchesne E (2008). The intrinsic functional organization of the
  brain is altered in autism. *NeuroImage*, 39(4), 1877–1885.
- Pelphrey KA, Morris JP, McCarthy G (2005). Neural basis of eye gaze processing
  deficits in autism. *Brain*, 128(5), 1038–1048.
- Lombardo MV et al. (2011). Specialization of right temporo-parietal junction
  for mentalizing and its relation to social impairments in autism. *NeuroImage*,
  56(3), 1832–1838.
- Uddin LQ, Menon V (2009). The anterior insula in autism: under-connected and
  under-examined. *Neuroscience & Biobehavioral Reviews*, 33(8), 1198–1203.

---

### 2.7 Alzheimer's Disease (AD)

**Seeds:** Hippocampus (HO Subcortical), entorhinal cortex (HO Cortical
parahippocampal gyrus), PCC / precuneus (HO Cortical), angular gyrus (HO
Cortical), nucleus basalis of Meynert (HO Subcortical — approximate; basal
forebrain atlas recommended).
**Networks:** DMN, Memory network.

AD follows a stereotyped spatiotemporal progression of tau pathology beginning
in the entorhinal cortex and hippocampus (Braak & Braak, 1991), which are also
the structures showing earliest volumetric loss and connectivity disruption.
A critical observation linking resting-state fMRI to AD pathophysiology is the
spatial congruence between the DMN and the pattern of amyloid-β deposition
(Buckner et al., 2005) — regions with the highest metabolic activity at rest
accumulate amyloid preferentially. Greicius et al. (2004) demonstrated reduced
DMN connectivity in clinically diagnosed AD, and subsequent work has shown DMN
disruption in preclinical and MCI stages (Sorg et al., 2007). The PCC/precuneus
seed for the DMN is particularly important in AD given that PCC hypometabolism
(on FDG-PET) is among the earliest functional biomarkers. The nucleus basalis
of Meynert, the primary source of cortical cholinergic innervation, shows early
neuronal loss in AD and its connectivity with cortical DMN nodes is a candidate
biomarker of cholinergic integrity.

**References:**
- Braak H, Braak E (1991). Neuropathological stageing of Alzheimer-related
  changes. *Acta Neuropathologica*, 82(4), 239–259.
- Buckner RL et al. (2005). Molecular, structural, and functional characterization
  of Alzheimer's disease: evidence for a relationship between default activity,
  amyloid, and memory. *Journal of Neuroscience*, 25(34), 7709–7717.
- Greicius MD, Srivastava G, Reiss AL, Menon V (2004). Default-mode network
  activity distinguishes Alzheimer's disease from healthy aging: evidence from
  functional MRI. *PNAS*, 101(13), 4637–4642.
- Sorg C et al. (2007). Selective changes of resting-state networks in individuals
  at risk for Alzheimer's disease. *PNAS*, 104(47), 18760–18765.
- Mesulam MM, Mufson EJ, Levey AI, Wainer BH (1983). Cholinergic innervation of
  cortex by the basal forebrain: cytochemistry and cortical connections of the
  septal area, diagonal band nuclei, nucleus basalis, and hypothalamus in the
  rhesus monkey. *Journal of Comparative Neurology*, 214(2), 170–197.

---

### 2.8 Frontotemporal Dementia (FTD)

FTD is not a single disease but a group of syndromes whose network involvement
varies systematically by clinical variant. The seminal observation of Seeley et
al. (2009) — that distinct neurodegenerative diseases target distinct large-scale
brain networks — forms the conceptual foundation for the network-based approach
taken here.

#### 2.8a Behavioural variant FTD (bvFTD)

**Seeds:** Anterior insula (HO Cortical), dACC (HO Cortical), frontal pole,
OFC (HO Cortical), right IFG.
**Networks:** Salience, Frontoparietal.

bvFTD preferentially degrades the Salience Network, beginning in the anterior
insula and dACC (Seeley et al., 2009; Seeley et al., 2007). This explains the
core clinical features: loss of empathy, disinhibition, and compulsive
behaviours reflecting failure of salience detection and social-emotional
processing.

#### 2.8b Semantic variant PPA (svPPA)

**Seeds:** Anterior temporal lobe / temporal pole (HO Cortical), amygdala
(HO Subcortical), parahippocampal gyrus (HO Cortical).
**Networks:** Semantic/ventral temporal, Limbic.

svPPA shows anterior temporal atrophy with semantic memory breakdown.
The temporal pole and anterior fusiform are primary targets; amygdala
involvement accounts for the prosopagnosia and socio-emotional deficits
frequently co-occurring.

#### 2.8c Non-fluent/agrammatic variant PPA (nfvPPA)

**Seeds:** Left IFG / Broca's area (BIP Language Atlas + HO Cortical),
premotor cortex, SMA (HO Cortical superior frontal).
**Networks:** Language, SMN.

nfvPPA targets the left dorsal language stream and premotor speech areas,
producing effortful, agrammatic speech with relative preservation of
comprehension.

#### 2.8d Logopenic variant PPA (lvPPA)

**Seeds:** Left posterior superior temporal gyrus / Wernicke's area
(HO Cortical), left inferior parietal lobule / angular gyrus (HO Cortical).
**Networks:** Language (posterior), DMN.

lvPPA involves the posterior language network and IPL, with phonological working
memory deficits as the cardinal feature. It is frequently associated with
underlying AD pathology.

**References (FTD section):**
- Seeley WW et al. (2009). Neurodegenerative diseases target large-scale human
  brain networks. *Neuron*, 62(1), 42–52.
- Seeley WW et al. (2007). Dissociable intrinsic connectivity networks for
  salience processing and executive control. *Journal of Neuroscience*, 27(9),
  2349–2356.
- Gorno-Tempini ML et al. (2011). Classification of primary progressive aphasia
  and its variants. *Neurology*, 76(11), 1006–1014.
- Rascovsky K et al. (2011). Sensitivity of revised diagnostic criteria for the
  behavioural variant of frontotemporal dementia. *Brain*, 134(9), 2456–2477.

---

### 2.9 Parkinson's Disease (PD)

**Seeds:** Putamen — caudal-motor subdivision (FSL Striatum), thalamus —
primary motor subdivision (FSL Thalamus), SMA (HO Cortical superior frontal),
precentral gyrus / M1 (HO Cortical), cerebellum motor lobules (SUIT atlas —
to be added), subthalamic nucleus (FSL STN atlas — to be added).
**Networks:** SMN, Frontostriatal (motor loop), Cerebellar-cortical.

PD is characterised by loss of dopaminergic neurons in the substantia nigra
pars compacta (SNc), leading to depletion of striatal dopamine and disruption
of the cortico-striato-thalamo-cortical motor loop. The posterior putamen
(sensorimotor striatum) is the first and most severely affected striatal region
(Kish et al., 1988), explaining the motor predominance of early PD. The FSL
Striatum caudal-motor subdivision captures this territory. Resting-state studies
demonstrate reduced putamen–motor cortex connectivity in PD, correlating with
motor symptom severity (Helmich et al., 2010). The subthalamic nucleus (STN) —
the target of DBS in advanced PD — is a critical node in the indirect and
hyperdirect pathways and shows pathological synchronisation with cortical motor
areas in PD (Brown et al., 2001); its inclusion requires a dedicated STN atlas
(e.g., the DISTAL atlas or FSL STN) not currently in the pipeline registry.
The cerebellum shows compensatory hyperactivity in early PD and its
connectivity with motor cortex via the thalamus is a target of interest for
distinguishing PD from atypical parkinsonian syndromes (Wu & Hallett, 2013).

**References:**
- Kish SJ, Shannak K, Hornykiewicz O (1988). Uneven pattern of dopamine loss in
  the striatum of patients with idiopathic Parkinson's disease. *New England
  Journal of Medicine*, 318(14), 876–880.
- Helmich RC, Derikx LC, Bakker M, Scheeringa R, Bloem BR, Toni I (2010).
  Spatial remapping of cortico-striatal connectivity in Parkinson's disease.
  *Cerebral Cortex*, 20(5), 1175–1186.
- Brown P et al. (2001). Dopamine dependency of oscillations between subthalamic
  nucleus and pallidum in Parkinson's disease. *Journal of Neuroscience*, 21(3),
  1033–1038.
- Wu T, Hallett M (2013). The cerebellum in Parkinson's disease. *Brain*, 136(3),
  696–709.

---

### 2.10 Hippocampal Sclerosis (HS)

**Seeds:** Hippocampus — ipsilateral and contralateral (HO Subcortical;
hippocampal subfield atlas recommended: ASHS or FreeSurfer 7 subfields),
parahippocampal gyrus / entorhinal cortex (HO Cortical), amygdala (HO
Subcortical), PCC (HO Cortical, DMN anchor).
**Networks:** MTL memory network, DMN, Limbic.

Hippocampal sclerosis — defined by selective loss of CA1 pyramidal neurons with
astrogliosis — is the most common pathological substrate of drug-resistant mesial
temporal lobe epilepsy (MTLE). Beyond the structural lesion, HS disrupts the
broader MTL memory network and propagates functional abnormalities to remote
regions via epileptic networks. Resting-state fMRI reveals reduced ipsilateral
hippocampal connectivity with the DMN and paradoxical increases in connectivity
reflecting ictal or interictal discharge propagation (Liao et al., 2010; Bettus
et al., 2009). Lateralisation of hippocampal connectivity asymmetry has
prognostic value for post-surgical memory outcome (Bonilha et al., 2010), making
asymmetric seed analysis (left vs right separately) essential in this profile.
Hippocampal subfields (CA1, CA3, CA4, dentate gyrus, subiculum) have distinct
vulnerability profiles in HS — CA1 is most affected — and subfield-level
connectivity analysis requires the ASHS atlas or FreeSurfer's hippocampal
subfields module, which are not yet integrated in the current registry.

**References:**
- Liao W et al. (2010). Default mode network abnormalities in mesial temporal
  lobe epilepsy: a study combining fMRI and DTI. *Journal of Neurology,
  Neurosurgery and Psychiatry*, 81(3), 328–336.
- Bettus G et al. (2009). Decreased basal fMRI functional connectivity in
  epileptogenic networks and contralateral compensatory mechanisms. *Human Brain
  Mapping*, 30(5), 1580–1591.
- Bonilha L et al. (2010). Presurgical connectome and postsurgical seizure
  control in temporal lobe epilepsy. *Neurology*, 81(19), 1704–1710.

---

## 3. Atlas Registry Summary

| Atlas | Source | Structures | Notes |
|---|---|---|---|
| Harvard-Oxford Cortical | FSL / nilearn | All major cortical gyri | 25% probability threshold, 2mm |
| Harvard-Oxford Subcortical | FSL / nilearn | Hippocampus, amygdala, caudate, putamen, NAcc, thalamus | 25% probability threshold, 2mm |
| FSL Striatum (7 subdivisions) | FSL | Limbic, executive, rostral-motor, caudal-motor, parietal, occipital, temporal striatum | Connectivity-based, 2mm |
| FSL Thalamus (7 subdivisions) | FSL | Primary motor, sensory, occipital, pre-frontal, pre-motor, posterior parietal, temporal | Connectivity-based, 2mm |
| BIP Language Atlas | Local | Surgical-grade language territories | Left hemisphere |
| Schaefer 2018 (400P, 17N) | nilearn | All cortical networks (for RSN masks) | Functional atlas — RSN masks only |
| Smith 2009 RSN | nilearn | 10 canonical RSNs | Functional atlas — RSN masks only |
| SUIT Cerebellar | To be added | Cerebellar lobules | Required for PD, ADHD |
| ASHS / FS Hippo Subfields | To be added | CA1–4, DG, subiculum | Required for HS |
| STN / DISTAL | To be added | Subthalamic nucleus | Required for PD |
| Basal Forebrain / nbM | To be added | Nucleus basalis of Meynert | Required for AD (cholinergic) |

---

## 4. General Methodological Notes

**SBA seed specificity:** Seeds are intentionally kept small and anatomically
focal. Larger seeds average across functionally heterogeneous subregions, washing
out the FC pattern of interest. Multiple seeds per network are run independently;
their SBA maps form the basis for the union mask passed to within-network ICA.

**Multiple atlases per structure:** Where a single atlas does not adequately
capture a structure (e.g., striatal subdivisions not parcellated in HO), seeds
from multiple atlases are unioned prior to timeseries extraction.

**Subject-specific refinement:** Atlas seeds are applied in MNI152NLin2009cAsym
space, to which each subject's BOLD has been normalised by fMRIPrep. For
subcortical structures with high inter-individual anatomical variability (NAcc,
hippocampus, amygdala), mri_synthseg-derived subject-specific parcellations
augment the atlas mask for RSN mask construction (not for SBA seeds, where
fixed anatomical definitions are preferred for cross-subject comparability).

**Lateralisation:** For profiles involving lateralised pathology (HS, language
tumours, nfvPPA), seeds are run separately for each hemisphere and results
reported with explicit lateralisation indices.

**N = 10 normative database:** All z-score comparisons are made against the
10-subject healthy volunteer normative sample. With n = 10, voxelwise SD
estimates are unreliable. ROI-level summary z-scores (mean Fisher-Z within
anatomical subregions) should be prioritised in clinical reports, and all
findings should be interpreted with explicit acknowledgement of this limitation.
