from reaction_mechanizer.pathway.reaction import *
from reaction_mechanizer.drawing.mechanism_reaction_visualizer import *

reac = ReactionMechanism.str_to_mechanism(
"""S+E->C
C->E+P""")

reac.set_rate_constants([{"kf":0.1,"kr":0},{"kf":0.05}])

vis = SimpleStepVisualizer(reac)
vis.progress_reaction({"S": 1, "E": 2,"C":0,"P":0}, 1000, 5000,events=[(200, ReactionEvent.CHANGE_CONCENTRATION, ("S", 0))],out = "out.png")