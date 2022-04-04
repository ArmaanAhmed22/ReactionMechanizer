![logo](readme_assets/reaction_mechanizer_logo.png)

---
# Reaction Mechanizer
Reaction Mechanizer is a Python tool that can be used to simulate chemical reactions.
# Theory
The rates at which steps in a reaction occur are governed by differential equations. For example, consider the equilibrium step:

```A+B<--->C+D```

Without loss of generality, if we focus solely on species `C`, notice that the rate of production of `C` is jointly proportional to the concentrations of `A` and `B` (more likely for the forward reaction to go through as the amount of interactions increase). However, the rate is also dictated by the reverse reaction, which is proportional to the interactions between `C` and `D`. So, the overall rate law for `C` is:

```dC/dt=k1*A*B - k2*C*D```

A system of ODEs can be built for any step, and in fact any full mechanism (a combination of steps).
# Documentation
To do
# Requirements
Reaction Mechanizer requires Python 3.10 in addition to the following libraries:
- matplotlib
- seaborn
- scipy
- pandas
- numpy

Specific version requirements are in `requirements.txt`
# Installation
This version is a pre-release, so Reaction Mechanizer can only be manually installed from source.

```
git clone https://github.com/ArmaanAhmed22/ReactionMechanizer.git & cd ReactionMechanizer & pip install .
```
# Usage
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

From an initial state, the step/mechanism progression can be seen using the `StepVisualizer` class as follows:

```python
from reaction_mechanism.drawing.mechanism_reaction_visualizer import StepVisualizer, ReactionEvent
vis = StepVisualizer(rm)  #or StepVisualizer(ss)
dataframe = vis.progress_reaction({"A": 1, "X": 0, "Y":0, "Z": 0, "B":0.1}, 1000, 5000, events=[(200, ReactionEvent.CHANGE_CONCENTRATION, ("A", 1))], out = "out.png")
"""Arguments:
1) initial state of species
2) time
3) granularity (how many times to evaluate the differential equations to derive our answer)
4) [optional] any events (ie adding certain species) that occur in the middle of the reaction (here, at time=200, with an increase of concentration of 1 for "A")
5) Where to save output file"""
```