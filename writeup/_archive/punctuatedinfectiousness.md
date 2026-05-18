# The impact of punctuated infectiousness on infectious disease dynamics 

## Introduction

The trajectory of an epidemic is governed by two key factors: the basic reproduction number, $R_0$, or the expected number of secondary infections that an infectious persion will produce in an otherwise susceptible population; and the generation interval distribution, $g(\tau)$, or the distribution of times at which those secondary infections are expected to occur, relative to the index case's infection ($\tau = 0$). 

Both the reproduction number and the generation interval distribution are expected values, averaged over a large population. At the individual level, the individual reproduction number [[x](https://www.nature.com/articles/nature04153)] is the expected number of new infections that a single infectious person will produce, such that $R_0 = \sum_i \nu_i$. Similarly, the generation interval distribution is the average of the times a single person is expected to transmit disease, such that $g(\tau) = \sum_i a_i(\tau)$, where $a_i(\tau)$ is the "individual infectiousness kernel" [[x](https://pubmed.ncbi.nlm.nih.gov/26674948/), [x](https://www.sciencedirect.com/science/article/pii/S0025556406002094)]. 

For many pathogens, the shape of the individual infectiousness kernel is poorly understood. Contact tracing to count infections is hard enough; measuring precisely the times of those infections is even harder. When it is done, it is usually averaged over many individuals, so that the "infectiousness period" assigned to a person in models reflects some average over multiple individuals -- either the generation interval itself, or some intermediate between the individual infectiousness kernel and the generation interval distribution. Some evidence suggests that the infectious period for SARS-CoV-2 and influenza, at the individual level, may be highly punctuated, such that opportuntiies for superspreading only exist over the course of a few hours. 

At the individual level, communicating about the likely timing of a person's infectiousness is important: it's good to know if you'll be very infectious for a short time, or mildly infectious for a long time, so that you can shape your behavior accordingly. Additionally, a punctuated individual infectiousness kernel could impact epidemic dynamics. Even if the epidemic has the same mean-field behavior, it could imapct the stochastic dynamics via highly synchronized transmission events, especially when an outbreak is small, It could interplay with control measures, especially those designed to intervene on inectiousness at the individual level like test- or symptom-based screening and isolation. Indeed, the impact of the punctuatedness of the individual infectiousness kernel has been examined as a sensitivity anlaysis in various epidemiological studies. (Dan; others?). However, we lack a comprehensive understanding of how the sharpenss of the individual infectiousness kernel, independent of other features like the individual reproduction number, basic reproduction number, and generation interval distribution, impacts epidemic dynamics under both uncontrolled and controlled settings. Furthermore, we lack an understanding of how to measure the shape of the individual infectiousness kernel reliably, along with the settings when the kernel's shape is statistically identifiable. 

Here, we address these gaps. 

## Expressing the Hudson-Kermack-McKendrick model in terms of individual infectiousness 

### The individual infectiousness profile 

Early work by Hudson, and later by Kermack and McKendrick [[x](https://link.springer.com/article/10.1007/BF02464423)], modeled disease transmission in a population using an integral equation. In modern notation [[x](https://www.tandfonline.com/doi/pdf/10.1080/17513758.2012.716454)], the force of infection is 

<a id="eq_foi"></a>
```math
F(t) = \int_0^\infty F(t-\tau) S(t-\tau) A(\tau) d\tau \tag{1}
```

where $S(t)$ is the density of susceptible individuals in the population at time $t$, and $A(\tau)$ is the population-level "infectiousness profile", a curve describing the expected contribution to the force of infection from an individual who was infected $\tau$ time units ago. The incidence of disease is simply the product of the force of infection and the susceptible density, $F(t) \cdot S(t)$; thus, the infectiousness profile $A(\tau)$ is thus the fundamental object that determines how an epidemic unfolds. 

The infectiousness profile is related to the basic reproduction number, $R_0$, and the generation interval, $g(\tau)$ [[x](https://www.tandfonline.com/doi/pdf/10.1080/17513758.2012.716454)]: 

$$ R_0 = \int_0^\infty A(\tau) d\tau $$ 

$$ g(\tau) = \frac{A(\tau)}{R_0} $$ 

Furthermore, particular choices of $A(\tau)$ yield disease transmission dynamics that are equivalent to more familiar ordinary differential equation-based models [[x](https://www.tandfonline.com/doi/pdf/10.1080/17513758.2012.716454)]. For example, For example, the SEIR model 

$$
\begin{align}
&\frac{dS}{dt} = -\beta S I  \\ 
&\frac{dE}{dt} = \beta S I - \gamma E  \\ 
&\frac{dI}{dt} = \gamma E - \alpha I  \\ 
&\frac{dR}{dt} = \alpha I  
\end{align}
$$

is obtained when 

$$ A(\tau) = \beta \frac{\gamma}{\gamma - \alpha} (e^{-\alpha \tau} - e^{-\gamma \tau}) $$

An individual person's infectiousness kernel $a_i(\tau)$ may differ substantially from $A(\tau)$, and multiple different $a_i(\tau)$ may yield the same $A(\tau)$. For example, the SEIR model is often conceptually motivated by considering an individual who undergoes a latent period of length $X \sim \text{Exp}(\gamma)$ during which they are not infectious, followed by an infectious period of length $Y \sim \text{Exp}(\alpha)$ during which they are infectious at level $\beta$. Mathematically, 

$$ a_i(\tau) = \begin{cases} 
\beta & \qquad \tau \in [\tau_i^{\text{on}}, \tau_i^{\text{on}} + \tau_i^{\text{off}}] \\
0 & \qquad \text{otherwise }
\end{cases} \qquad \text{ where } \tau_i^{\text{on}} \sim \text{Exp}(\gamma) \text{ and } \tau_i^{\text{off}} \sim \text{Exp}(\alpha)$$ 

It can be shown that $E_i [a_i(\tau)] = \lim_{n \rightarrow \infty} \frac{1}{n} \sum_{i = 1}^n a_i(\tau) = A(\tau)$ [[x](https://www.tandfonline.com/doi/pdf/10.1080/17513758.2012.716454)]. However, other $a_i(\tau)$ also yield the same $A_i(\tau)$, including 

$$ a_i(\tau) = \beta \frac{\gamma}{\gamma - \alpha} (e^{-\alpha \tau} - e^{-\gamma \tau}) $$ 

where each person's infectiousness profile is a smooth curve that is identical to the population-level infectiousness profile, and 

$$ a_i(\tau) = \frac{\beta}{\gamma} \cdot \delta_{t_i}(\tau) \qquad \text{ where } t_i \text{ is distributed according to the generation interval density, } A(\tau)/R_0$$ 

where each person's infectiousness profile is a delta function, such that all their infectiousness is concentrated at a single moment (**Figure XX**). 

Like $A(\tau)$, $a_i(\tau)$ is related to the individual reproduction number, $\nu_i$, and the individual generation interval distribution, $\xi_i(\tau)$: 

$$ \nu_i = \int_\tau=0^\infty a_i(\tau) $$ 

$$ \xi_i(\tau) = a_i(\tau) / \nu_i $$ 


### A Gamma convolutional model for the individual infectiousness profile 

The "stepwise" expression for $a_i(\tau)$ in Eq XX differs from the "smooth" and "spike" expressions in Eqs XX-XX in a fundamental way: For Eq XX, $\nu_i$ differs across individuals ($\nu_i \sim  XX$), while for Eq XX - XX, $\nu_i$ is uniform across individuals ($\nu_i = R_0$). We will briefly compare how the unconstrained disease transmission dynamics differ for Eq XX - XX, but because of this fundamental dissimilarity -- and because the impact of variation in $\nu_i$ on epidemic dynamics has been discussed exstensively elsewhere [x] -- we will focus most of the manuscript on individual infectiousness profiles of the type in Eq XX - XX, where the individual infectiousness $\nu_i$ is uniform across the population. That way, the impact of the timing of infectiousness can be examined in isolation. 

To achieve this, we seek a family of functions that allow us to model individual infectiousness profiles with the following characteristics: 

* All have uniform total infectiousness: $\int_0^\infty a_i(\tau) d\tau = R_0$ 
* Smooth, single-peaked distributions 
* All have identical shape, except for a time shift (no mixtures of broad smooth distributions with narrow spiky distributions)

A candidate set of functions that obeys this is the set of Gamma convolutional curves. Let $A(\tau)$ be Gamma distributed (up to a scaling constant): 

$$
A(\tau) = R_0 f_A(\tau) \qquad \text{ where } f_A(\tau) \sim \text{Gamma}(\alpha, r), \qquad r = \alpha / \mu
$$

where $\mu$ is the mean generation interval and $\alpha$ controls the shape of the population-level infectiousness profile. 

Now, let 

$$
a_i(\tau) = R_0 f_a(\tau - s_i) \qquad \text{ where } f_a(\tau) \sim \text{Gamma}(\kappa, r) \qquad \text{ and } s_i \sim \text{Gamma}(\alpha - \kappa, r)
$$

Here, $\kappa \in (0, \alpha)$ is a parameter governing the punctuatedness of $a_i(\tau)$, with smaller $\kappa$ yielding a more concentrated infectiousness profile (**Figure XX**). Since the sum of two Gamma-distributed random variables with the same rate is also Gamma distributed (with the same rate and with shape equal to the sum of the two component shapes), this formulation is guaranteed to converge to $A(\tau)$ in expectation (**Supplementary Methods**). Furthermore, the $k$-values may differ across individuals; regardless of the distribution they are drawn from, the Gamma sum property ensures that $E_i[a_i(\tau)] = A(\tau)$ (**Supplementary Methods**, **Figure XX**).

### Splitting the infectiousness profile into biological infectiousness and contacts 

From equations XX-XX, the individual infectiousness profile can be expressed in terms of the individual generation interval distribution and the individual reproduction number: 

$$ a_i(\tau) = \xi_i(\tau) \nu_i $$ 

A different and more general decomposition is also possible, into factors relating to a person's "biological infectiousness", $b_i(\tau)$, and "contact potential", $c_i(t)$: 

$$ a_i(\tau) = b_i(\tau) \, c_i(t_i + \tau) $$

The biological infectiousness, $b_i(\tau)$, is a nonnegative density function that integrates to 1 ($\int_0^\infty b_i(\tau) d\tau = 1$) that describes how a person's potential infectiousness is spread over time. The contact potential, $c_i(t)$, is a nonnegative time-varying function that describes the expected number of individuals person $i$ would infect if all their infectiousness were concentrated at the single moment $t$ (i.e., if $b_i(\tau)$ were a delta function). Thus, $c_i(t)$ incorporates both a person's contact rate and the probability of transmission given contact. It can immediately be seen that if $c_i(t)$ is constant, then we recover Eq XX where $b_i(\tau)$ is the individual generation interval distribution and $c_i(t) = \nu_i$. For time-varying $c_i(t)$, this form shows that the individual reproduction number and the individual generation interval distribution depends not just on the time since infection $\tau$, but also on clock time $t$. For infection time $t^*$, we obtain **CHECK THESE**

$$ \nu_i(t^*) = \int_{\tau = 0}^{\infty} b_i(\tau) \, c_i(t^* + \tau) \, d\tau $$

$$ \xi_i(t^*, \tau) = \frac{b_i(\tau) \, c_i(t^* + \tau)}{\nu_i(t^*)} $$ 

This decomposition feeds back into the population-level model. If we assume that $b_i$ and $c_i$ are independent (*i.e.*, an individual's biological infectiousness is independent of their contact pattern), then 

$$
A_i(\tau) = E[a_i(\tau)] = B(\tau) \cdot C(t_i + \tau)
$$

where $B(\tau) = E_i[b_i(\tau)]$ and $C(t+i + \tau) = E_i[c_i(t + \tau)]$. Likewise, if C is contstant in time, this reduces to the familiar 

$$
A_i(\tau) = g(\tau) R_0
$$

Splitting $a_i(\tau)$ into biological infectiousness $b_i$ and contact potential $c_i$ allows us to separately consider how the shape of a person's infectiousness, unconditional on their contact patterns, impacts transmission dynamics under various contact patterns. 

For the gamma convolutional model, we have 

$$
b_i(\tau) = f_a(\tau - s_i)
$$

and 

$$
c_i(t) = R_0 \text{ for all } i,t
$$


## The impact of the individual infectiousness profile on uncontrolled epidemic dynamics 

### Constant contacts 

* Stepwise, smooth, and spike all yield the same mean dynamics 
* Stepwise has higher extinction probability, explained by overdispersion
* Spike has later and more variable establishment properties, including time-to-n-infected and peak, and exponential growth rate is more variable. We can show this through simualtions and through analytic approaches, following Morris' time shift ideas [[x](https://link.springer.com/article/10.1007/s00285-024-02132-6)]

With constant contacts, kappa doesn't affect extinction or overdispersion — only growth timing

### Time-varying contacts

but introduce time-varying contacts and suddenly kappa controls overdispersion

For the rest of our analysis, we focus on individual infectiousness profiles that can be expressed using the Gamma convolutional framework so that each person's total infectiousness is uniform across the population and the only thing that varies is the punctuation and timing of the individual infectiousness profile.

* Spike emphasizes variation in contacts, leading to overdispersion. 

## The impact of punctuated infectiousness on outbeak interventions



### Detection and isolation 

* Detection and isolation is more effective for punctuated kernels 
* Detect-and-isolate policies shift the overall generation interval distribution in a way that depends on $\kappa$ - leading to a partially self-defeating interventino. 

### Contact tracing 
* Amidst all-same and kappa-variable scenarios, with tracing delays. 

### Gathering size restrictions 

* More effective for smooth profiles, since they reduce mean transmission without losing the overdispersion that was suppressing outbreak probability for spiky profiles. 

## Identifiability of punctuated infectiousness 

While the concentration of the individual infectiousness profile can impact how an epidemic unfolds, we lack frameworks for estimating it. Many approaches exist for estimating the generation interval, but there has not been a systematic approach developed to assess how individual infectiousness varies and to specify where on the continuum we sit -- i.e., how much of the variation in the generation interval distribution is attributable to variation between individuals, and how much to variation in infection timing within individuals. 

# Discussion 

# References 

# Acknowledgements 

# Funding 

# Author contributions 

# Supplementary Information 

## Supplementary Methods


## Supplementary Figures 


## Supplementary Tables


<!--
## Notes on writeup structure and uncontrolled dynamics section (Claude review, 2026-02-27)

## Structure

The overall flow is clean and logical: framework → uncontrolled dynamics → interventions → identifiability. A few observations:

1. **The bio/contact decomposition** (Section "Splitting the infectiousness profile...") is introduced in the framework before it's needed. It's essential for the overdispersion story (sampling vs. averaging), which lives in the uncontrolled dynamics section. You could either (a) keep it as-is (all tools upfront) or (b) introduce it at the point where time-varying contacts enter the picture. Either is defensible, but (b) would give the reader a clearer motivation for why the decomposition matters.

2. **The uncontrolled dynamics section** would benefit from two subsections that mirror the decomposition:
   - **Constant contacts** (Findings 1–7): same mean-field dynamics, stepwise vs. smooth/spike extinction, min-of-order-statistics growth delay, burstiness, within-individual correlation structure
   - **Time-varying contacts** (Findings 10–13): the sampling-vs-averaging mechanism, overdispersion from punctuated × periodic contacts, epidemic consequences

   This makes the punchline of the first subsection ("with constant contacts, kappa doesn't affect extinction or overdispersion — only growth timing") a natural setup for the second subsection ("but introduce time-varying contacts and suddenly kappa controls overdispersion").

3. **The interventions** have a nice progression from simpler to more complex. The unified detect-and-isolate framework (independent ↔ anchored triggers) followed by contact tracing (two-generation coupling) followed by gathering restrictions (variance mechanism) reads well.

## Additions for uncontrolled dynamics

Your current bullet points cover: mean dynamics, extinction probability, establishment timing/growth rate, and overdispersion from contacts. The user asked specifically about **peaks**. Here's my assessment:

### Peaks — yes, worth adding

**Peak timing distribution.** Since the establishment delay for spiky profiles is ~1 day (at R0=2) and grows with R0, this should shift peak timing later by a similar amount. More importantly, the *variance* in peak timing should be higher for spiky profiles (inheriting the higher variance in establishment timing from Finding 4). But once established, mean-field dynamics dominate, so **peak height should converge across kappa** — the difference is when, not how high.

This would be a compact, satisfying result: "The epidemic peak occurs at approximately the same height regardless of kappa, but the timing of the peak is more variable for punctuated profiles, reflecting the delayed and more stochastic establishment phase."

It's essentially a downstream consequence of Finding 4 (establishment timing variance), so it doesn't require much new machinery — just extending the epidemic simulations to track peak incidence timing and height.

### Other potential additions

**Epidemic duration / tail behavior.** Does the epidemic tail (the last 10% of infections) differ across kappa? In the smooth case, the tail should be more predictable; in the spike case, late stragglers whose onset shifts fall late could create more variable epidemic tails. This is minor but would round out the picture.

**Serial interval ≠ generation interval.** If you later connect to observation models, it's worth noting that the *observed* serial interval distribution could differ from the generation interval even when A(tau) is the same, because serial intervals depend on the relationship between infectiousness timing and symptom timing. With biologically anchored symptoms (the rho=1 scenario from Section 14), spiky profiles would produce more clustered serial intervals (all transmissions from one infector have nearly identical serial intervals), while smooth profiles would produce more dispersed serial intervals. This is more of an identifiability-relevant observation, so it might fit better in the identifiability section.

### What I'd skip

**Network/spatial effects** — these are important in general but open a large can of worms that would dilute the paper's focus on the temporal structure of infectiousness.

**Phylogenetic signatures** — Finding 7 (within-individual correlation structure) hints at this, but developing it fully would require a phylodynamic simulation framework that seems out of scope.

### Summary recommendation for the uncontrolled section

The four bullets you have, plus **peak timing/height** and the **two-subsection split** (constant contacts → time-varying contacts), would make this section complete. The peak result is the most natural addition — it's what a reader would expect to see and it falls out directly from existing simulations with minimal new analysis.
-->
