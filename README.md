# Size ratios in continuous chain model
Continuous chain model allows to reseive an analytical expresions for size characteristics like gyration radius and hydrodynamic radius. The down side is that the expressions are rather big and thus hard to evaluate by had when taken from a published paper. This project allows to get the numbers from those expressions. 

## Parameters that can be calculated
1. Squared gyration radius $\langle R_g^2\rangle = \frac{1}{2N}\sum_{i=1}^N\sum_{j=1}^N \langle(\vec{r}_i-\vec{r}_j)^2\rangle$ where $\vec{r}$ are the radius vectors of monomers and $\langle\ldots\rangle$ defines ansamble averaging. In general it can be written in a form $\langle R_g^2 \rangle = A(\mathcal{F})l^2N^{2\nu}$ where $A(\mathcal{F})$ is a universal amplitude that depends on polymer architecture but does not depend on the chemistry or detail of the model, $\nu$ is a Flory exponent and $l$ is a Khun length that depend on the chemistry or detail of the model. This app outputs a gyration radius $\sqrt{\langle R_g^2\rangle}$ assuming $l=1$ for 
- Ideal or gaussian polymer chain
- Polymer in good solution with excluded volume interaction (real chain). The calculation are done using [Douglas-Freed approximation](https://pubs.acs.org/doi/abs/10.1021/ma00141a026).

2. Hydrodynamic radius  $\langle R_h^{-1}\rangle = \frac{1}{N}\sum_{i=1}^N\sum_{j=1}^N \langle(\vec{r}_i-\vec{r}_j)^2{-1}\rangle$ for an ideal chain. Again it can be in general presented as $\langle R_h \rangle = A(\mathcal{F})l N^{\nu}$ with $l$ in the output assumed as 1.

3. Size ratio $\rho=\frac{ \langle R_g^2 \rangle _{branched}}{\langle R_g^2\rangle _{linear}}$ that compairs the sizes of molecules with different architectures of the same mass. The output again is provided for:
-ideal chain
-real chain

4. Size ratio $\rho=\frac{R_g}{R_h}$ for ideal chains. This ratio compairs different measurements for the same molecular architecture.

## Architectures avaliable

1. Rosette polymers: $f_c$ open chains and $f_r$ closed chains (rings) connected to a single branching point. Results used for this project where published in:
- Gyration radius for an ideal chain and corresponding $g_c$ ratio [publication](https://iopscience.iop.org/article/10.1088/1751-8113/48/13/135001), [ArXiv](https://arxiv.org/abs/1412.2553)
- Hydrodinamic for an ideal chain and corresponding $\rho$ ratio [publication](https://www.nature.com/articles/s41598-020-70649-z) avaliable in open access
- Gyration radius for a real chain and corresponding $g_c$ ratio[publication](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.105.034502)
2. Pom-pom polymers: a backbone chain with $f_1$ and $f_2$ open chains attached to each of its ends. Results used for this project where published in:
- Gyration radius and corresponding $g_c$ ratio: Symmetric $f_1=f_2$ [publication](https://www.sciencedirect.com/science/article/abs/pii/S0167732221001823?via%3Dihub), [ArXiv](https://arxiv.org/abs/2012.00469), Assymetric [publication open acess](https://icmp.lviv.ua/journal/zbirnyk.110/23302/abstract.html), influence of the backbone [publication open acess](https://icmp.lviv.ua/journal/zbirnyk.114/23301/abstract.html)
- Hydrodynamic radius and corresponding $\rho$ ratio [publication](https://iopscience.iop.org/article/10.1088/1751-8121/ac5508), [ArXiv](https://arxiv.org/abs/2201.09053)
3. Dumbbell polymers: a backbone chain with a closed chain attached to each of its ends
- Gyration radius and corresponding $g_c$ ratio [publication](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.108.034502), [ArXiv](https://arxiv.org/abs/2303.11062)
- Hydrodynamic radius and corresponding $\rho$ ratio [publication](https://iopscience.iop.org/article/10.1088/1751-8121/ac5508), [ArXiv](https://arxiv.org/abs/2201.09053)
4. Snowflake polymers: $f_s$ chains are connected in a branching point with $f-1$ open chains attached to the oposite end of eah of then
- Gyration radius and corresponding $g_c$ ratio [publication](https://www.sciencedirect.com/science/article/pii/S0167732223022365?via%3Dihub) avaliable in open access
- Hydrodynamic radius expressions were derived for this app. For the purposeds of citatios please cite [publication](https://iopscience.iop.org/article/10.1088/1751-8121/ac5508) as it was used as a base.
5. Bottlebrush polymers: a long backbone with a number of branching points evenly spaced along it. Each branching point contains $f$ open chains
- Ideal chains [publication](https://iopscience.iop.org/article/10.1088/1751-8121/ac5508), [ArXiv](https://arxiv.org/abs/2201.09053)
- Real chains -- publication is in the process of preparation 

## Files and How to use
This project contains:
- a folder with schematics of the architectures considered. All figures were made for this project and do not apper in the papers
- a single python file that contains both the programed expressions as well as interface. It requres tkinter and numpy to run 
- a executable file for windows that you can run without instalation

**Note** for both .py and .exe a folder with pictures next to the file is required.   

With any questions regarding the analytical calculations you can contact me on [ResearchGate](https://www.researchgate.net/profile/Khristine-Haydukivska)