from typing import Dict, List
import pytest
from reaction_mechanizer.pathway.reaction import ReactionMechanism, SimpleStep
from reaction_mechanizer.drawing.mechanism_reaction_visualizer import ReactionVisualizer

epsilon = 0.001
end_time = 1000
granularity = 1000


@pytest.mark.parametrize("k_values, eq_step, initial_conditions, expected_K", [
    ((kf := 1, kr := 1), SimpleStep.str_to_step("A + B -> C"), [{"A": 100, "B": 200, "C": 0}, {"A": 100, "B": 200, "C": 300}], kf/kr),
    ((kf := 2, kr := 0.25), SimpleStep.str_to_step("A + B -> C"), [{"A": 100, "B": 2, "C": 0}, {"A": 100, "B": 200, "C": 3000}], kf/kr)
])
def test_equilibrium_step(k_values: tuple, eq_step: SimpleStep, initial_conditions: List[Dict[str, float]], expected_K: float):
    eq_step.set_rate_constant(kf=k_values[0], kr=k_values[1])
    vis = ReactionVisualizer(eq_step)
    K_error_list = []
    for initial_condition in initial_conditions:
        df = vis.progress_reaction(initial_condition, end_time, granularity)
        ends = df.iloc[-1]
        numerator = 1
        denominator = 1
        for key, i in eq_step.products.items():
            numerator *= ends[key]**i
        for key, i in eq_step.reactants.items():
            denominator *= ends[key]**i
        K_error_list.append(abs(numerator/denominator - expected_K))
    for K_error in K_error_list:
        if K_error > epsilon:
            assert False
    assert True


@pytest.mark.parametrize("k_values_list, mechanism, initial_condition, species_of_importance, expected_concentration", [
    (
        [{"kf": 1, "kr": 0.05}, {"kf": 0.2}],
        ReactionMechanism.str_to_mechanism("""S+E->C
                                           C->E+P"""),
        {"S": 2, "E": 1, "C": 0, "P": 0},
        "S",
        0
     ),
    (
        [{"kf": 1, "kr": 0.05}, {"kf": 0.2}],
        ReactionMechanism.str_to_mechanism("""S+E->C
                                           C->E+P"""),
        {"S": 2, "E": 1, "C": 0, "P": 0},
        "P",
        2
     ),
    (
        [{"kf": 1, "kr": 0.05}, {"kf": 0.2}],
        ReactionMechanism.str_to_mechanism("""S+E->C
                                           C->E+P"""),
        {"S": 2, "E": 1, "C": 0.5, "P": 0},
        "P",
        2.5
     ),
])
def test_mechanism(
        k_values_list: List[Dict[str, float]],
        mechanism: ReactionMechanism,
        initial_condition: Dict[str, float],
        species_of_importance: str,
        expected_concentration):
    mechanism.set_rate_constants(k_values_list)
    vis = ReactionVisualizer(mechanism)
    df = vis.progress_reaction(initial_condition, end_time, granularity)
    assert abs(df.iloc[-1][species_of_importance] - expected_concentration) <= epsilon


@pytest.mark.parametrize("k_values_list, mechanism, initial_condition, species_of_importance, expected_concentration", [
    (
        [{"kf": 1, "kr": 0.05}, {"kf": 0.2}],
        ReactionMechanism.str_to_mechanism("""S+E->C
                                           C->E+P"""),
        {"S": 2, "E": 1, "C": 0, "P": 0},
        "S",
        0
     ),
    (
        [{"kf": 1, "kr": 0.05}, {"kf": 0.2}],
        ReactionMechanism.str_to_mechanism("""S+E->C
                                           C->E+P"""),
        {"S": 2, "E": 1, "C": 0, "P": 0},
        "P",
        2
     ),
    (
        [{"kf": 1, "kr": 0.05}, {"kf": 0.2}],
        ReactionMechanism.str_to_mechanism("""S+E->C
                                           C->E+P"""),
        {"S": 2, "E": 1, "C": 0.5, "P": 0},
        "P",
        2.5
     ),
])
def test_mechanism_robust(
        k_values_list: List[Dict[str, float]],
        mechanism: ReactionMechanism,
        initial_condition: Dict[str, float],
        species_of_importance: str,
        expected_concentration):
    mechanism.set_rate_constants(k_values_list)
    vis = ReactionVisualizer(mechanism)
    df = vis.progress_reaction_robust(initial_condition, end_time, method="Radau")
    assert abs(df.iloc[-1][species_of_importance] - expected_concentration) <= epsilon
