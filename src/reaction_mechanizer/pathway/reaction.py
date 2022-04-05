from types import LambdaType
from typing import Any, Dict, Union
import re
import typing


class SimpleStep:

    def __init__(self, reactants: Dict[str, float], products: Dict[str, float]):
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
        if "->" in str_step:
            str_step = str_step.replace("\t", "").replace(" ", "")
            reac_str, prod_str = str_step.split("->")
            return SimpleStep.str_to_step(reac_str, _optional_flag=True) + SimpleStep.str_to_step(prod_str, _optional_flag=False)

        if "+" in str_step:
            return typing.cast('SimpleStep', sum([SimpleStep.str_to_step(component, _optional_flag=_optional_flag) for component in str_step.split("+")]))

        match = re.match(r"^[0-9]+\/[0-9]+|^[0-9]+", str_step)
        cutoff: int = match.span()[1] if match is not None else (1, str_step := "1"+str_step)[0]

        if "/" in str_step:
            parts: list[str] = str_step[:cutoff].split("/")
            coef = int(parts[0]) / int(parts[1])
        else:
            coef = int(str_step[:cutoff])

        if _optional_flag is True:
            return SimpleStep({str_step[cutoff:]: coef}, {})
        else:
            return SimpleStep({}, {str_step[cutoff:]: coef})

    def set_rate_constant(self, kf: float = 0, kr: float = 0) -> None:
        self.kf = kf
        self.kr = kr

    def get_differential_equations(self) -> dict[str, 'DifferentialEquationModel']:
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

    def get_differential_equation_of(self, thing: str) -> 'DifferentialEquationModel':
        arguments: str = ",".join(self.reactants.keys())
        products = [thing for thing in self.products.keys() if thing not in self.reactants.keys()]
        if len(products) != 0:
            arguments += ","+",".join(products)
        reactant_multiplication: str = "*".join([f"""kwargs["{thing}"]**{coef}""" for thing, coef in self.reactants.items()])
        product_multiplication: str = "*".join([f"""kwargs["{thing}"]**{coef}""" for thing, coef in self.products.items()])

        ode = DifferentialEquationModel("")
        if thing in self.reactants:
            coef_reactant = self.reactants[thing]
            ode = DifferentialEquationModel.sum_differential_equations(
                [ode, DifferentialEquationModel(f"-{coef_reactant}*{self.kf}*{reactant_multiplication} + {coef_reactant}*{self.kr}*{product_multiplication}")])
        if thing in self.products:
            coef_product = self.products[thing]
            new_ode = DifferentialEquationModel(f'{coef_product}*{self.kf}*{reactant_multiplication} - {self.kr}*{coef_product}*{product_multiplication}')
            ode = DifferentialEquationModel.sum_differential_equations([ode, new_ode])
        return ode

    def set_rate_constant_from_K(self, K, normalizing_kr=1) -> None:
        self.kr = normalizing_kr
        self.kf = self.kr * K

    def __add__(self, other: 'SimpleStep') -> 'SimpleStep':
        def dict_addition(dict_from, add_to): return {**dict_from, **{thing: dict_from.get(thing, 0)+add_to[thing] for thing in add_to.keys()}}
        return SimpleStep(dict_addition(self.reactants, other.reactants), dict_addition(self.products, other.products))

    def __radd__(self, other: Union['SimpleStep', int]) -> 'SimpleStep':
        if type(other) is int:
            return self
        return self.__add__(typing.cast('SimpleStep', other))

    def __str__(self):
        reactants: str = " + ".join([str(coef) + thing for thing, coef in self.reactants.items()])
        products: str = " + ".join([str(coef) + thing for thing, coef in self.products.items()])
        return reactants + " -> " + products

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.reactants == other.reactants and self.products == other.products
        else:
            return False


class ReactionMechanism:

    class MechanismException(Exception):
        pass

    def __init__(self, steps: list[SimpleStep], strictness: str = "warning"):
        """_summary_

        Args:
            steps (list[SimpleStep]): _description_
            strictness (str, optional): Levels: "warning", "error". Defaults to "warning".
        """
        self.steps = steps

    @staticmethod
    def str_to_mechanism(str_mechanism: str, /, strictness: str = "warning"):
        steps = [SimpleStep.str_to_step(str_step) for str_step in str_mechanism.split("\n")]
        return ReactionMechanism(steps)

    def set_rate_constants(self, rates: list[dict[str, float]]):
        for new_rate, step in zip(rates, self.steps):
            step.set_rate_constant(**new_rate)

    def get_differential_equations(self) -> dict[str, 'DifferentialEquationModel']:
        out_ode_prev: dict[str, list['DifferentialEquationModel']] = {}
        out_ode: dict[str, 'DifferentialEquationModel'] = {}
        for step in self.steps:
            for key, ode in step.get_differential_equations().items():
                if key in out_ode_prev:
                    out_ode_prev[key].append(ode)
                else:
                    out_ode_prev[key] = [ode]

        for key in out_ode_prev.keys():
            out_ode[key] = DifferentialEquationModel.sum_differential_equations(out_ode_prev[key])
        return out_ode

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
        """Create a first order Differential Equation representation

        Args:
            model_str (str): Should involve all the variables using "kwargs['var']".
        """
        self.model_str = model_str

    def get_lambda(self) -> LambdaType:
        return typing.cast(LambdaType, eval(f"lambda **kwargs: {self.model_str}"))

    @staticmethod
    def sum_differential_equations(list_ode: list['DifferentialEquationModel']) -> 'DifferentialEquationModel':
        return DifferentialEquationModel('+'.join([ode.model_str for ode in list_ode]))
