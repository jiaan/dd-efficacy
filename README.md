# Efficacy of dynamical decoupling #
(Working title for now; can we find a better one?) 

Joint work by Jiaan Qi, Xu Xiansong, Dario Poletti, and Hui Khoon Ng

## Goal of the work ## 
To understand the efficacy of dynamical decoupling in removing noise in quantum devices, using fault tolerance ideas.
First, a bit of background on fault tolerance. Standard fault-tolerant quantum computing makes use of quantum error correction (QEC) codes to provide resilience against noise. It is the only known route to large-scale, useful quantum computers. In particular, fault tolerance theory discusses a few things: The use of QEC = increased complexity in the quantum computation = more things to go wrong. How to accomplish QEC when the QEC operations themselves are noisy -- how to remove noise, rather than add noise.A protocol for scaling up the noise protection by putting in more physical resources into the QEC (to increase the error correction ability of the QEC code), for increased computational accuracy.A threshold - known as the “quantum accuracy threshold” - for the noise below which the protocol works (i.e., can scale up the resource use to increase accuracy).Break-even point (or “pseudothresholds”) for the noise for which a FIXED SCALE protection works (i.e., better than no QEC). The quantum accuracy threshold can (sort of) be thought of as the minimum over the break-even points for all scales.
Draw parallels with this to analyze dynamical decoupling (DD). Replace “QEC” above by “DD”. For DD, scaling up = CDD.

## Main points of our work ##
Break-even points: Presence of a break-even point (for stated figure of merit: error phase or infidelity) for single level DD (PDD as our example).Analytical calculationNumerical estimation for more general noise modelsIllustrate the sensitivity of the conclusions on the choice of figure of meritNo quantum accuracy threshold:As CDD scales up, the correction capability does not scale up enough to overcome the increase in complexity - no accuracy thresholdThis holds even for ideal CDD, i.e., without noisy pulses, due to the increased correction power not being sufficient to overcome the increase in time taken and hence the accumulated error Why CDD fails? Partly (the other part is just accumulated leftover higher-order errors) due to the growing effective bath noise — the alpha, which leads to faster and fast noise dynamics that the DD pulse rate eventually fail to keep up with. Understanding this point would be very interesting (though to check how much of this is already in the Lidar paper).For noisy CDD, this is even worse. 

##Main messages
Presence of a break-even point for the use of dynamical decoupling at each level of CDD. After some level of concatenation, break-even point no longer nonzero - no use going further in the concatenation.Our work shows how to estimate the optimal CDD level for a given experimental situation. (Although we don’t know \alpha??)Why no quantum accuracy threshold? We show why CDD fails eventually even for the ideal case, and identify features of the DD scaling-up scheme one must have to overcome the accumulated error -- hints at future directions to investigate.

## Possible referee questions:
What about Uhrig decoupling? What is the break-even point for that?Other schemes of scaling up?
