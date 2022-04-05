from enum import Enum
from typing import Any, Tuple, Union
import typing
from reaction_mechanizer.pathway.reaction import DifferentialEquationModel, ReactionMechanism, SimpleStep
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


class ReactionEvent(Enum):
    CHANGE_CONCENTRATION = 0
    SET_CONCENTRATION = 1
    SMOOTH_CHANGE_CONCENTRATION = 2


class StepVisualizer:

    def __init__(self, simple_step: Union[SimpleStep, ReactionMechanism]):
        self.simple_step: Union[SimpleStep, ReactionMechanism] = simple_step

    def get_states(self,
                   initial_state: dict[str, float],
                   time_end: float,
                   number_steps: int,
                   initial_time: float = 0,
                   ode_override: Union[dict[str, DifferentialEquationModel], None] = None) -> np.ndarray[Any, Any]:
        ode_dict: dict[str, DifferentialEquationModel] = self.simple_step.get_differential_equations()
        ode_dict_temp = {key: ode_dict[key] for key in initial_state.keys()}
        ode_dict = ode_dict_temp
        if ode_override is not None:
            for key, ode in ode_override.items():
                ode_dict[key] = ode
        ode_function = _get_simple_step_ode_function(ode_dict, list(initial_state.keys()))

        times = np.linspace(initial_time, time_end, number_steps)

        cur_state = list(initial_state.values())

        return typing.cast(np.ndarray[Any, Any], odeint(ode_function, cur_state, times))

    def progress_reaction(self,
                          initial_state: dict[str, float],
                          time_end: float, number_steps: int,
                          events: Union[list[Tuple[float, ReactionEvent, Tuple[Any]]], None] = None,
                          out: Union[str, None] = None):
        data: np.ndarray[Any, Any] = np.ndarray((0, 0))
        if events is not None:
            sorted_events = (*sorted(events, key=lambda x: x[0]), (time_end, None, tuple()))
            cur_state = dict(initial_state)
            prev_time_point_discretized: float = 0

            for time_point, reaction_event_type, additional_info in sorted_events:
                # First, discretize the time points so that everything is measured to the granularity of "number_steps"
                time_interval = time_end / number_steps
                time_point_discretized = round(time_point / time_interval) * time_interval

                cur_number_steps = round((time_point - prev_time_point_discretized) / time_end * number_steps)

                cur_data: np.ndarray[Any, Any] = self.get_states(cur_state, time_point, cur_number_steps, initial_time=prev_time_point_discretized)
                data = cur_data if len(data) == 0 else typing.cast(np.ndarray[Any, Any], np.concatenate([data, cur_data]))

                if reaction_event_type == ReactionEvent.CHANGE_CONCENTRATION:
                    for index, key in enumerate(cur_state.keys()):
                        if key == additional_info[0]:
                            cur_state[key] = cur_data[-1, index] + additional_info[1]
                        else:
                            cur_state[key] = cur_data[-1, index]
                elif reaction_event_type == ReactionEvent.SET_CONCENTRATION:
                    for index, key in enumerate(cur_state.keys()):
                        if key == additional_info[0]:
                            cur_state[key] = cur_data[-1, index] + additional_info[1]
                        else:
                            cur_state[key] = cur_data[-1, index]
                elif reaction_event_type == ReactionEvent.SMOOTH_CHANGE_CONCENTRATION:
                    pass
                prev_time_point_discretized = time_point_discretized
        else:
            data = self.get_states(initial_state, time_end, number_steps)
        times = np.linspace(0, time_end, number_steps)

        fig, ax = plt.subplots()
        plt.tight_layout()
        if out:
            for i, thing in enumerate(initial_state.keys()):
                sns.lineplot(x=times, y=data[:, i], label="$"+thing+"$", ax=ax)
            ax.legend()
            sns.despine(ax=ax)
            ax.margins(x=0, y=0)
            _, top = ax.get_ylim()
            ax.set_ylim([0, top*1.05])
            plt.savefig(str(out), bbox_inches="tight", dpi=600)
        return pd.DataFrame({thing: data[:, i] for i, thing in enumerate(initial_state.keys())})


def _get_simple_step_ode_function(differential_equations: dict[str, DifferentialEquationModel], state_order: list[str]):
    active_odes = {key: val.get_lambda() for key, val in differential_equations.items()}

    def simple_step_ode_function(cur_state, t):
        cur_state_order = state_order
        dict_state = {thing: val for thing, val in zip(cur_state_order, cur_state)}
        out = []
        for _, ode in active_odes.items():
            out.append(ode(**dict_state))
        return out
    return simple_step_ode_function
