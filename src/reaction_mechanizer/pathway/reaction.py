from collections import Counter
from types import LambdaType
from typing import Any, Dict, Union
import re


class SimpleStep:

    def __init__(self, reactants: Union[Dict[str, float], Counter], products: Union[Dict[str, float], Counter]):
        self.liquid_solid_reactants: Counter = Counter()
        self.reactants: Counter = Counter()
        for thing, coef in reactants.items():
            if "(l)" in thing or "(s)" in thing:
                self.liquid_solid_reactants[thing] = coef
            else:
                self.reactants[thing] = coef
        self.products: Counter = Counter(products)
        self.liquid_solid_products: Counter = Counter()

        for thing, coef in products.items():
            if "(l)" in thing or "(s)" in thing:
                self.liquid_solid_products[thing] = coef
            else:
                self.products[thing] = coef

        self.kf: float = 0
        self.kr: float = 0

    @staticmethod
    def str_to_step(str_step: str, _optional_flag: Any = None) -> 'SimpleStep':
        if "->" in str_step:
            reac_str, prod_str = str_step.split("->")
            return SimpleStep.str_to_step(reac_str, _optional_flag=True) + SimpleStep.str_to_step(prod_str, _optional_flag=False)

        if "+" in str_step:
            return sum([SimpleStep.str_to_step(component, _optional_flag=_optional_flag) for component in str_step.replace(" ", "").split("+")])

        match = re.match(r"^[0-9]+\/[0-9]+|^[0-9]+", str_step)
        cutoff: int = match.span()[1] if match is not None else (1, str_step := "1"+str_step)[0]

        if "/" in str_step:
            parts: list[str, str] = str_step[:cutoff].split("/")
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

    def get_differential_equations(self) -> dict[str, LambdaType]:
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
            print(f"d{thing}/dt = {raw_func}")
            out[thing] = (eval(raw_func), raw_func)
        return out_ode

    def get_differential_equation_of(self, thing: str) -> dict[str, LambdaType]:
        arguments: str = ",".join(self.reactants.keys())
        products = [thing for thing in self.products.keys() if thing not in self.reactants.keys()]
        if len(products) != 0:
            arguments += ","+",".join(products)
        reactant_multiplication: str = "*".join([f"""kwargs["{thing}"]**{coef}""" for thing, coef in self.reactants.items()])
        product_multiplication: str = "*".join([f"""kwargs["{thing}"]**{coef}""" for thing, coef in self.products.items()])

        out_raw_str = ""
        if thing in self.reactants:
            coef_reactant = self.reactants[thing]
            out_raw_str += f"lambda **kwargs: -{coef_reactant}*{self.kf}*{reactant_multiplication} + {coef_reactant}*{self.kr}*{product_multiplication}"
        if thing in self.products:
            coef_product = self.products[thing]
            new_lambda_str = f'lambda **kwargs: {coef_product}*{self.kf}*{reactant_multiplication} - {self.kr}*{coef_product}*{product_multiplication}'
            out_raw_str = f"lambda **kwargs: ({out_raw_str})(**kwargs) + ({new_lambda_str})(**kwargs)" if not out_raw_str == "" else new_lambda_str
        return eval(out_raw_str)

    def set_rate_constant_from_K(self, K, normalizing_kr=1) -> None:
        self.kr = normalizing_kr
        self.kf = self.kr * K

    def __add__(self, other: 'SimpleStep') -> 'SimpleStep':
        return SimpleStep(self.reactants + other.reactants, self.products + other.products)

    def __radd__(self, other: Union['SimpleStep', int]) -> 'SimpleStep':
        if type(other) is int:
            return self
        return self.__add__(other)

    def __str__(self):
        reactants: str = " + ".join([str(coef) + thing for thing, coef in self.reactants.items()])
        products: str = " + ".join([str(coef) + thing for thing, coef in self.products.items()])
        return reactants + " -> " + products


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
        """if not (kind_mistakes:=self.mistakes(),)[0] is None:
            match strictness:
                case "warning":
                    print("WARNING: mismatch of the following\n"+str(kind_mistakes))
                case "error":
                    raise self.MechanismException("ERROR: mismatch of the following\n"+str(kind_mistakes))"""

    @staticmethod
    def str_to_mechanism(str_mechanism: str, /, strictness: str = "warning"):
        steps = [SimpleStep.str_to_step(str_step) for str_step in str_mechanism.split("\n")]
        return ReactionMechanism(steps)

    def set_rate_constants(self, rates: list[dict]):
        for new_rate, step in zip(rates, self.steps):
            step.set_rate_constant(**new_rate)

    def _get_differential_equations(self):
        out: dict[str, LambdaType] = {}
        for step in self.steps:
            for key, ode in step.get_differential_equations().items():
                if key in out:
                    print("YES")
                    prev_ode = out[key]

                    out[key] = lambda **kwargs: prev_ode(**kwargs) + ode(**kwargs)
                else:
                    out[key] = ode
        return out

    def get_differential_equations(self):
        out_ode_prev = {}
        out_ode = {}
        for step in self.steps:
            for key, ode in step.get_differential_equations().items():
                if key in out_ode_prev:
                    out_ode_prev[key].append(ode)
                else:
                    out_ode_prev[key] = [ode]

        for key in out_ode_prev.keys():
            out_ode[key] = DifferentialEquationModel.sum_differential_equations(out_ode_prev[key])
            print(out_ode[key].model_str)
        return out_ode

    def mistakes(self) -> Counter | None:

        epsilon = 0.00001

        cur_mistakes = Counter()
        for step in self.steps:
            for thing, coef in step.reactants:
                cur_mistakes[thing] = cur_mistakes.get(thing, 0) + coef
            for thing, coef in step.products:
                cur_mistakes[thing] = cur_mistakes.get(thing, 0) - coef

        for thing, error in cur_mistakes.items():
            if abs(error) < epsilon:
                cur_mistakes.pop(thing)
        return cur_mistakes if len(cur_mistakes) >= 1 else None

    def __str__(self) -> str:
        return "\n".join([str(step) for step in self.steps])


class DifferentialEquationModel:
    def __init__(self, model_str: str):
        """Create a first order Differential Equation representation

        Args:
            model_str (str): Should involve all the variables using "kwargs['var']".
        """
        self.model_str = model_str

    def get_lambda(self) -> LambdaType:
        return eval(f"lambda **kwargs: {self.model_str}")

    @staticmethod
    def sum_differential_equations(list_ode: list['DifferentialEquationModel']) -> 'DifferentialEquationModel':
        print(list_ode)
        return DifferentialEquationModel('+'.join([ode.model_str for ode in list_ode]))
