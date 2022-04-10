![logo](readme_assets/reaction_mechanizer_logo.png)

---
# Reaction Mechanizer
![Tests](https://github.com/ArmaanAhmed22/ReactionMechanizer/actions/workflows/tests.yml/badge.svg)

Reaction Mechanizer is a Python tool that can be used to simulate chemical reactions.

# Table of Contents
1. [Theory](#theory)
2. [Documentation](#documentation)
3. [Requirements](#requirements)
4. [Installation](#installation)
5. [Usage](#usage)
# Theory <a id="theory"></a>
The rates at which steps in a reaction occur are governed by differential equations. For example, consider the equilibrium step:

```A+B<--->C+D```

Without loss of generality, if we focus solely on species `C`, notice that the rate of production of `C` is jointly proportional to the concentrations of `A` and `B` (more likely for the forward reaction to go through as the amount of interactions increase). However, the rate is also dictated by the reverse reaction, which is proportional to the interactions between `C` and `D`. So, the overall rate law for `C` is:

```dC/dt=k1*A*B - k2*C*D```

A system of ODEs can be built for any step, and any full mechanism (a combination of steps).
# Documentation <a id="documentation"></a>
Compiled documentation for this whole package is hosted [here](https://armaanahmed22.github.io/ReactionMechanizer/).
# Requirements <a id="requirements"></a>
Reaction Mechanizer requires Python 3.8 in addition to the following libraries:
- matplotlib
- seaborn
- scipy
- pandas
- numpy

Specific version requirements are in `requirements.txt`
# Installation <a id="installation"></a>
This version is a pre-release, so Reaction Mechanizer can only be manually installed from source.

```bash
git clone https://github.com/ArmaanAhmed22/ReactionMechanizer.git && cd ReactionMechanizer && pip install .
```
# Usage <a id="usage"></a>
Reaction Mechanizer can either simulate a `SimpleStep` or a `ReactionMechanism`. Either can be created from a string:

```python
from reaction_mechanism.pathway.reaction import SimpleStep, ReactionMechanism

ss = SimpleStep.str_to_step("2A+1/2B->1/3C+D")

rm = ReactionMechanism.str_to_mechanism(
"""A->X+Y
X->Z
Z+Y->B"""
)
```

From an initial state, the step/mechanism progression can be seen using the `ReactionVisualizer` class as follows:

```python
from reaction_mechanism.drawing.mechanism_reaction_visualizer import ReactionVisualizer, ReactionEvent
vis = ReactionVisualizer(rm)  #or ReactionVisualizer(ss)
dataframe = vis.progress_reaction({"A": 1, "X": 0, "Y":0, "Z": 0, "B":0.1}, 1000, 5000, events=[(200, ReactionEvent.CHANGE_CONCENTRATION, ("A", 1))], out = "out.png")
"""Arguments:
1) initial state of species
2) time
3) granularity (how many times to evaluate the differential equations to derive our answer)
4) [optional] any events (ie adding certain species) that occur in the middle of the reaction (here, at time=200, with an increase of concentration of 1 for "A")
5) Where to save output file"""
```