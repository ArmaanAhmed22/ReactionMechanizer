from abc import ABC, abstractmethod
from enum import Enum
from math import ceil
from typing import Tuple
from reaction_mechanizer.pathway.reaction import ReactionMechanism, SimpleStep
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

    def __init__(self, simple_step: SimpleStep|ReactionMechanism):
        self.simple_step: SimpleStep|ReactionMechanism = simple_step

    def get_states(self, initial_state: dict, time_end, number_steps, initial_time: int = 0, ode_override: dict = None):
        ode_dict = self.simple_step.get_differential_equations()
        ode_dict_temp = {key: ode_dict[key] for key in initial_state.keys()}
        ode_dict = ode_dict_temp
        if ode_override is not None:
            for key, ode in ode_override.items():
                ode_dict[key] = ode
        ode_function = _get_simple_step_ode_function(ode_dict, initial_state.keys())

        times = np.linspace(initial_time, time_end, number_steps)

        cur_state = list(initial_state.values())

        return odeint(ode_function, cur_state, times)

    def progress_reaction(self, initial_state: dict, time_end, number_steps, events: list[Tuple[float, ReactionEvent, Tuple]] = None, out: str = None):
        if events is not None:
            sorted_events = (*sorted(events, key=lambda x: x[0]), (time_end, None, tuple()))
            cur_state = dict(initial_state)
            prev_time_point_discretized = 0

            data = None

            for time_point, reaction_event_type, additional_info in sorted_events:
                # First, discretize the time points so that everything is measured to the granularity of "number_steps"
                time_interval = time_end / number_steps
                time_point_discretized = round(time_point / time_interval) * time_interval

                cur_number_steps = round((time_point - prev_time_point_discretized) / time_end * number_steps)

                cur_data = self.get_states(cur_state, time_point, cur_number_steps, initial_time=prev_time_point_discretized)
                data = cur_data if data is None else np.concatenate((data, cur_data))
                match reaction_event_type:
                    case ReactionEvent.CHANGE_CONCENTRATION:
                        for index, key in enumerate(cur_state.keys()):
                            if key == additional_info[0]:
                                cur_state[key] = cur_data[-1, index] + additional_info[1]
                            else:
                                cur_state[key] = cur_data[-1, index]
                    case ReactionEvent.SET_CONCENTRATION:
                        for index, key in enumerate(cur_state.keys()):
                            if key == additional_info[0]:
                                cur_state[key] = cur_data[-1, index] + additional_info[1]
                            else:
                                cur_state[key] = cur_data[-1, index]
                    case ReactionEvent.SMOOTH_CHANGE_CONCENTRATION:
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
            _,top = ax.get_ylim()
            ax.set_ylim([0,top*1.05])
            plt.savefig(str(out), bbox_inches="tight", dpi=600)
        return pd.DataFrame({thing: data[:, i] for i, thing in enumerate(initial_state.keys())})


def _get_simple_step_ode_function(differential_equations: dict, state_order: list):
    active_odes = {key: val.get_lambda() for key, val in differential_equations.items()}

    def simple_step_ode_function(cur_state, t):
        cur_state_order = state_order
        dict_state = {thing: val for thing, val in zip(cur_state_order, cur_state)}
        out = []
        for _, ode in active_odes.items():
            out.append(ode(**dict_state))
        return out
    return simple_step_ode_function
