from typing import Dict, List
from reaction_mechanizer.pathway.reaction import DifferentialEquationModel, SimpleStep
from reaction_mechanizer.drawing.mechanism_reaction_visualizer import ReactionMechanism
import pytest


@pytest.mark.parametrize("string_input, expected", [
    ("A+B->C+D", SimpleStep({"A": 1, "B": 1}, {"C": 1, "D": 1})),
    ("A +B -> C+      D", SimpleStep({"A": 1, "B": 1}, {"C": 1, "D": 1})),
    ("1/2A+3B->C+5/6D", SimpleStep({"A": 1/2, "B": 3}, {"C": 1, "D": 5/6})),
    ("1/2A+ 3B->C+  3/800D", SimpleStep({"A": 1/2, "B": 3}, {"C": 1, "D": 3/800}))
])
def test_str_to_step(string_input, expected):
    assert SimpleStep.str_to_step(string_input) == expected


@pytest.mark.parametrize("string_input, expected", [
    (
        """A->X
        X->C""", ReactionMechanism(
            [SimpleStep({"A": 1}, {"X": 1}),
             SimpleStep({"X": 1}, {"C": 1})]
        )
    ),
    (
        """A -> X
        X-> 1/5C""", ReactionMechanism(
            [SimpleStep({"A": 1}, {"X": 1}),
             SimpleStep({"X": 1}, {"C": 1/5})]
        )
    ),
    (
        """A+B+C->X
        X->C
        C->D+E""", ReactionMechanism(
            [SimpleStep({"A": 1, "B": 1, "C": 1}, {"X": 1}),
             SimpleStep({"X": 1}, {"C": 1}),
             SimpleStep({"C": 1}, {"D": 1, "E": 1})]
        )
    )
])
def test_str_to_mechanism(string_input, expected):
    assert ReactionMechanism.str_to_mechanism(string_input) == expected


@pytest.mark.parametrize("rate_of, k_values, coordinates, input_step, ode_expected_lambda", [
    ("A", (kf := 1, kr := 1), [{"A": 1, "B": 1, "C": 1, "D": 1}, {"A": 2, "B": 2, "C": 2, "D": 1}], "A+B->C+D", lambda A, B, C, D: -kf*A*B + kr*C*D),
    ("D", (kf := 1, kr := 1), [{"A": 1, "B": 1, "C": 1, "D": 1}, {"A": 2, "B": 2, "C": 2, "D": 1}], "A+B->C+D", lambda A, B, C, D: kf*A*B - kr*C*D)
])
def test_differential_equations_step(rate_of: str, k_values: tuple, coordinates: List[Dict[str, float]], input_step: str, ode_expected_lambda):
    matchting_rates = []
    step = SimpleStep.str_to_step(input_step)
    step.set_rate_constant(kf=k_values[0], kr=k_values[1])
    ode_input: DifferentialEquationModel = step.get_differential_equation_of(rate_of)
    for coord in coordinates:
        matchting_rates.append(ode_input.get_lambda()(**coord) == ode_expected_lambda(**coord))
    assert sum(matchting_rates) == len(coordinates)


@pytest.mark.parametrize("simple_step, expected", [
    (SimpleStep({"A": 1, "B": 2}, {"C": 1, "D": 3}), "1A + 2B -> 1C + 3D")
])
def test_str_simple_step(simple_step: SimpleStep, expected: str):
    assert str(simple_step) == expected


@pytest.mark.parametrize("reaction_mechanism, expected_intermediates", [
    (
        ReactionMechanism.str_to_mechanism(
            """A+B->2C
            C->D
            D+C->J"""),
        ["C", "D"]
    ),
    (
        ReactionMechanism.str_to_mechanism(
            """2A+B->C
            """),
        []
    )
])
def test_intermediates(reaction_mechanism: ReactionMechanism, expected_intermediates: List[str]):
    assert set(reaction_mechanism.get_intermediates()) == set(expected_intermediates)
