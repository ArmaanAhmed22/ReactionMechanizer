"""Contains tools to model a single step or a whole mechanism for a reaction.
"""
from types import LambdaType
from typing import Any, Dict, List, Union
import re
import typing


class SimpleStep:
    """SimpleStep represents a single step in a mechanism. For example, it can model aA+bB->cC+dD using a K constant or kf and kr rate constants
    """
    def __init__(self, reactants: Dict[str, float], products: Dict[str, float]):
        """Create a SimpleStep object

        Args:
            reactants (Dict[str, float]): List of reactants in the following form: {"A":a,"B":b}. To represent fractional coefficients, `a` and `b` are floats
            products (Dict[str, float]): List of products in the following form: {"A":a,"B":b}. To represent fractional coefficients, `a` and `b` are floats
        """
        self.liquid_solid_reactants: Dict[str, float] = {}
        self.reactants: Dict[str, float] = {}
        for thing, coef in reactants.items():
            if "(l)" in thing or "(s)" in thing:
                self.liquid_solid_reactants[thing] = coef
            else:
                if thing in self.reactants:
                    self.reactants[thing] += coef
                else:
                    self.reactants[thing] = coef
        self.products: Dict[str, float] = {}
        self.liquid_solid_products: Dict[str, float] = {}

        for thing, coef in products.items():
            if "(l)" in thing or "(s)" in thing:
                self.liquid_solid_products[thing] = coef
            else:
                if thing in self.products:
                    self.products[thing] += coef
                else:
                    self.products[thing] = coef

        self.kf: float = 0
        self.kr: float = 0

    @staticmethod
    def str_to_step(str_step: str, _optional_flag: Any = None) -> 'SimpleStep':
        """Turns a string representation of a reaction to a `SimpleStep` representation of that reaction.

        Args:
            str_step (str): \
                The string representation of the reaction of the form: "aA + bB -> cC + dD". `a`, `b`, `c`, and `d` are coefficents \
                    (if one is fractional, use fraction notation rather than decimal notation ["1/2" GOOD, "0.5" BAD]). \
                        "->" separates reactants and products, and "A", "B", "C", "D" are names of the species involved in the reaction \
                            (they shouldn't contain any special characters)
            _optional_flag (Any, optional): Dummy flag used internally to transfer information during recursive steps.

        Returns:
            SimpleStep: SimpleStep representation of `str_step`. `self.kf` and `self.kr` are set to 0
        """
        if "->" in str_step:
            str_step = str_step.replace("\t", "").replace(" ", "")
            reac_str, prod_str = str_step.split("->")
            return SimpleStep.str_to_step(reac_str, _optional_flag=True) + SimpleStep.str_to_step(prod_str, _optional_flag=False)

        if "+" in str_step:
            return typing.cast('SimpleStep', sum([SimpleStep.str_to_step(component, _optional_flag=_optional_flag) for component in str_step.split("+")]))

        match = re.match(r"^[0-9]+\/[0-9]+|^[0-9]+", str_step)
        cutoff: int = match.span()[1] if match is not None else (1, str_step := "1"+str_step)[0]

        if "/" in str_step:
            parts: List[str] = str_step[:cutoff].split("/")
            coef = int(parts[0]) / int(parts[1])
        else:
            coef = int(str_step[:cutoff])

        if _optional_flag is True:
            return SimpleStep({str_step[cutoff:]: coef}, {})
        else:
            return SimpleStep({}, {str_step[cutoff:]: coef})

    def set_rate_constant(self, kf: float = 0, kr: float = 0) -> None:
        """Set the rate constants for `self.kf` and `self.kr`

        Args:
            kf (float, optional): The kf. Defaults to 0.
            kr (float, optional): The kr. Defaults to 0.
        """
        self.kf = kf
        self.kr = kr

    def get_differential_equations(self) -> Dict[str, 'DifferentialEquationModel']:
        """Get system of first order differential equations model for the progression of this step.

        Returns:
            Dict[str, DifferentialEquationModel]:
                Get dictionary consisting of all the species in this step as keys and all the corresponding differential equation models as values.
        """
        arguments: str = ",".join(self.reactants.keys())
        products = [thing for thing in self.products.keys() if thing not in self.reactants.keys()]
        if len(products) != 0:
            arguments += ","+",".join(products)
        reactant_multiplication: str = "*".join([f"""kwargs["{thing}"]**{coef}""" for thing, coef in self.reactants.items()])
        product_multiplication: str = "*".join([f"""kwargs["{thing}"]**{coef}""" for thing, coef in self.products.items()])

        out = {}
        out_raw_str = {}
        out_ode = {}

        for thing, coef in self.reactants.items():
            out_raw_str[thing] = f"lambda **kwargs: -{coef}*{self.kf}*{reactant_multiplication} + {coef}*{self.kr}*{product_multiplication}"
            out_ode[thing] = DifferentialEquationModel(f"-{coef}*{self.kf}*{reactant_multiplication} + {coef}*{self.kr}*{product_multiplication}")
        for thing, coef in self.products.items():
            if thing in out_raw_str:
                new_lambda_str = f'lambda **kwargs: {coef}*{self.kf}*{reactant_multiplication} - {self.kr}*{coef}*{product_multiplication}'
                out_raw_str[thing] = f"lambda **kwargs: ({out_raw_str[thing]})(**kwargs) + ({new_lambda_str})(**kwargs)"
                new_ode = DifferentialEquationModel(f"{coef}*{self.kf}*{reactant_multiplication} - {self.kr}*{coef}*{product_multiplication}")
                out_ode[thing] = DifferentialEquationModel.sum_differential_equations([out_ode[thing], new_ode])
            else:
                out_raw_str[thing] = f"lambda **kwargs: {coef}*{self.kf}*{reactant_multiplication} - {self.kr}*{coef}*{product_multiplication}"
                out_ode[thing] = DifferentialEquationModel(f"{coef}*{self.kf}*{reactant_multiplication} - {self.kr}*{coef}*{product_multiplication}")
        for thing, raw_func in out_raw_str.items():
            out[thing] = (eval(raw_func), raw_func)
        return out_ode

    def get_differential_equation_of(self, species: str) -> 'DifferentialEquationModel':
        """Get a differential equation for just a single `species` in this step.

        Args:
            species (str): The species to get the differential equation of

        Returns:
            DifferentialEquationModel: The differential equation model of the `species`
        """
        arguments: str = ",".join(self.reactants.keys())
        products = [species for species in self.products.keys() if species not in self.reactants.keys()]
        if len(products) != 0:
            arguments += ","+",".join(products)
        reactant_multiplication: str = "*".join([f"""kwargs["{species}"]**{coef}""" for species, coef in self.reactants.items()])
        product_multiplication: str = "*".join([f"""kwargs["{species}"]**{coef}""" for species, coef in self.products.items()])

        ode = DifferentialEquationModel("")
        if species in self.reactants:
            coef_reactant = self.reactants[species]
            ode = DifferentialEquationModel.sum_differential_equations(
                [ode, DifferentialEquationModel(f"-{coef_reactant}*{self.kf}*{reactant_multiplication} + {coef_reactant}*{self.kr}*{product_multiplication}")])
        if species in self.products:
            coef_product = self.products[species]
            new_ode = DifferentialEquationModel(f'{coef_product}*{self.kf}*{reactant_multiplication} - {self.kr}*{coef_product}*{product_multiplication}')
            ode = DifferentialEquationModel.sum_differential_equations([ode, new_ode])
        return ode

    def set_rate_constant_from_K(self, K: float, normalizing_kr=1) -> None:
        """Set reaction rate using K-value normalized to a certain kr (default is 1). K is mathematically equivalent to kf/kr

        Args:
            K (float): The K-value
            normalizing_kr (int, optional): What kr to use. Defaults to 1.
        """
        self.kr = normalizing_kr
        self.kf = self.kr * K

    def __add__(self, other: 'SimpleStep') -> 'SimpleStep':
        def dict_addition(dict_from, add_to): return {**dict_from, **{thing: dict_from.get(thing, 0)+add_to[thing] for thing in add_to.keys()}}
        return SimpleStep(dict_addition(self.reactants, other.reactants), dict_addition(self.products, other.products))

    def __radd__(self, other: Union['SimpleStep', int]) -> 'SimpleStep':
        if type(other) is int:
            return self
        return self.__add__(typing.cast('SimpleStep', other))

    def __str__(self) -> str:
        """String representation of this `SimpleStep`. Uses the following format: "aA + bB -> cC + dD"

        Returns:
            str: string version
        """
        reactants: str = " + ".join([str(coef) + thing for thing, coef in self.reactants.items()])
        products: str = " + ".join([str(coef) + thing for thing, coef in self.products.items()])
        return reactants + " -> " + products

    def __eq__(self, other: Any) -> bool:
        """Check equality between two `SimpleStep`s. Equality is when the species and corresponding coefficients are the same among both objects.

        Args:
            other (Any): Other object.

        Returns:
            bool: whether the two objects are equal or not.
        """
        if isinstance(other, self.__class__):
            return self.reactants == other.reactants and self.products == other.products
        else:
            return False


class ReactionMechanism:
    """`ReactionMechanism` is a series of `SimpleStep`s that is analogous to a chemical mechanism behind a reaction.
    """

    class MechanismException(Exception):
        pass

    def __init__(self, steps: List[SimpleStep]):
        """Create a ReactionMechanism object

        Args:
            steps (List[SimpleStep]): list of simple steps to conjoin into a `ReactionMechanism`
        """
        self.steps = steps

    @staticmethod
    def str_to_mechanism(str_mechanism: str) -> 'ReactionMechanism':
        """Create `ReactionMechanism` from string representation of mechanism

        Args:
            str_mechanism (str): String mechanism. It is a list of `SimpleStep`s (following their string syntax), and each `SimpleStep` is separated by `\\n`

        Returns:
            ReactionMechanism: Resultant `ReactionMechanism`
        """
        steps = [SimpleStep.str_to_step(str_step) for str_step in str_mechanism.split("\n")]
        return ReactionMechanism(steps)

    def set_rate_constants(self, rates: List[Dict[str, float]]):
        """Set the `self.kf` and `self.kr` values for all the simple steps in this mechanism.

        Args:
            rates (List[Dict[str, float]]): Ordinal list of dictionaries containing entries for either kf, kr, both, or none.
        """
        for new_rate, step in zip(rates, self.steps):
            step.set_rate_constant(**new_rate)

    def get_differential_equations(self) -> Dict[str, 'DifferentialEquationModel']:
        """Get the set of `DifferentialEquationModel` representation of this mechanism

        Returns:
            Dict[str, DifferentialEquationModel]: Dictionary whose entries represent the differential equation model of each of the species in the mechanism.
        """
        out_ode_prev: Dict[str, List['DifferentialEquationModel']] = {}
        out_ode: Dict[str, 'DifferentialEquationModel'] = {}
        for step in self.steps:
            for key, ode in step.get_differential_equations().items():
                if key in out_ode_prev:
                    out_ode_prev[key].append(ode)
                else:
                    out_ode_prev[key] = [ode]

        for key in out_ode_prev.keys():
            out_ode[key] = DifferentialEquationModel.sum_differential_equations(out_ode_prev[key])
        return out_ode

    def get_intermediates(self):
        epsilon = 0.0001
        net_species: Dict[str, float] = {}
        for step in self.steps:
            for spec, coef in step.reactants.items():
                net_species[spec] = net_species.get(spec, 0) + coef
            for spec, coef in step.products.items():
                net_species[spec] = net_species.get(spec, 0) - coef
        return [spec for spec, number in net_species.items() if abs(number) < epsilon]

    def __str__(self) -> str:
        return "\n".join([str(step) for step in self.steps])

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if len(other.steps) == len(self.steps):
                return sum([other_step == cur_step for other_step, cur_step in zip(other.steps, self.steps)]) \
                       == len(self.steps)
            return False
        return False


class DifferentialEquationModel:
    def __init__(self, model_str: str):
        """Create an object representing a first order differential equation. It may contain multiple dependent variables.

        Args:
            model_str (str): Should involve all the variables using "kwargs['var']".
        """
        self.model_str = model_str

    def get_lambda(self) -> LambdaType:
        """Get the lambda version of this `DifferentialEquationModel`. Use keyword arguments to specify variable values.

        Returns:
            LambdaType: Lambda version of this `DifferentialEquationModel`
        """
        return typing.cast(LambdaType, eval(f"lambda **kwargs: {self.model_str}"))

    @staticmethod
    def sum_differential_equations(list_ode: List['DifferentialEquationModel']) -> 'DifferentialEquationModel':
        """Sum together the differential equations in the list

        Args:
            list_ode (List[DifferentialEquationModel]): List of `DifferentialEquationModel`s

        Returns:
            DifferentialEquationModel: The resultant differential equation
        """
        return DifferentialEquationModel('+'.join([ode.model_str for ode in list_ode]))
