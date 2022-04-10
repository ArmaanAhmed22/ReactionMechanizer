"""Contains tools to visualize a reaction.
"""
from enum import Enum
from typing import Any, Dict, List, Tuple, Union
import typing
from reaction_mechanizer.pathway.reaction import DifferentialEquationModel, ReactionMechanism, SimpleStep
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


class ReactionEvent(Enum):
    """Enum for the possible reaction events, such as:\n
    CHANGE_CONCENTRATION: how much the concentration of the species should be changed by.
        Additional Info: Tuple(str: species of interest, float: change in concentration)
    SET_CONCENTRATION: what to set the concentration of the species to.
        Additional Info: Tuple(str: species of interest, float: new concentration)
    SMOOTH_CHANGE_CONCENTRATION: TBD
    """
    CHANGE_CONCENTRATION = 0
    SET_CONCENTRATION = 1
    SMOOTH_CHANGE_CONCENTRATION = 2


class ReactionVisualizer:
    """Visualier for either `SimpleStep` or `ReactionMechanism`
    """
    def __init__(self, reaction: Union[SimpleStep, ReactionMechanism]):
        """Create a `ReactionVisualizer` object to model the progression of a reaction

        Args:
            reaction (Union[SimpleStep, ReactionMechanism]): The reaction object to model
        """
        self.reaction: Union[SimpleStep, ReactionMechanism] = reaction

    def get_states(self,
                   initial_state: Dict[str, float],
                   time_end: float,
                   number_steps: int,
                   initial_time: float = 0,
                   ode_override: Union[Dict[str, DifferentialEquationModel], None] = None) -> Any:
        """Get concentration of the species in this reaction, with model specifications given.

        Args:
            initial_state (Dict[str, float]): initial concentration of all species in reaction
            time_end (float): The end time for this model
            number_steps (int): The granularity of this model. The higher the number of steps, the more accurate the model.
            initial_time (float, optional): The time to start the model at. Defaults to 0.
            ode_override (Union[Dict[str, DifferentialEquationModel], None], optional): \
                Dictionary containing the species to override the differential equation of using the provided one. Defaults to None.

        Returns:
            Any: 2D array where the rows represent the concentrations of the species at different times \
                (between `initial_time` and `end_time` and using `number_steps`). The columns are the species in the order given by `initial_state`
        """
        ode_dict: Dict[str, DifferentialEquationModel] = self.reaction.get_differential_equations()
        ode_dict_temp = {key: ode_dict[key] for key in initial_state.keys()}
        ode_dict = ode_dict_temp
        if ode_override is not None:
            for key, ode in ode_override.items():
                ode_dict[key] = ode
        ode_function = _get_simple_step_ode_function(ode_dict, list(initial_state.keys()))

        times = np.linspace(initial_time, time_end, number_steps)

        cur_state = list(initial_state.values())

        return odeint(ode_function, cur_state, times)

    def progress_reaction(self,
                          initial_state: Dict[str, float],
                          time_end: float,
                          number_steps: int,
                          events: Union[List[Tuple[float, ReactionEvent, Tuple[Any]]], None] = None,
                          out: Union[str, None] = None) -> pd.DataFrame:
        """Generate model for reaction

        Args:
            initial_state (Dict[str, float]): initial concentration of all species in reaction
            time_end (float): The end time for this model
            number_steps (int): The granularity of this model. The higher the number of steps, the more accurate the model.
            events (Union[List[Tuple[float, ReactionEvent, Tuple[Any]]], None], optional): \
                The list of events to occur during a specified time in the reaction. \
                    A single event is represented by a tuple holding the time of the perturbation, the type of perturbation (`ReactionEvent`), \
                        and the additional information associated with the `ReactionEvent` selected. Defaults to None.
            out (Union[str, None], optional): \
                If a string is added, a png visually representing the reaction is created at the specified location (and a `DataFrame` is returned). \
                    Otherwise, just the `DataFrame` is returned. Defaults to None.

        Returns:
            pd.DataFrame: DataFrame representing the concentrations of the species in the reaction
        """
        data: Any = np.ndarray((0, 0))
        if events is not None:
            sorted_events = (*sorted(events, key=lambda x: x[0]), (time_end, None, tuple()))
            cur_state = dict(initial_state)
            prev_time_point_discretized: float = 0

            for time_point, reaction_event_type, additional_info in sorted_events:
                # First, discretize the time points so that everything is measured to the granularity of "number_steps"
                time_interval = time_end / number_steps
                time_point_discretized = round(time_point / time_interval) * time_interval

                cur_number_steps = round((time_point - prev_time_point_discretized) / time_end * number_steps)

                cur_data = self.get_states(cur_state, time_point, cur_number_steps, initial_time=prev_time_point_discretized)
                data = cur_data if len(data) == 0 else typing.cast(Any, np.concatenate([data, cur_data]))

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


def _get_simple_step_ode_function(differential_equations: Dict[str, DifferentialEquationModel], state_order: List[str]):
    active_odes = {key: val.get_lambda() for key, val in differential_equations.items()}

    def simple_step_ode_function(cur_state, t):
        cur_state_order = state_order
        dict_state = {thing: val for thing, val in zip(cur_state_order, cur_state)}
        out = []
        for _, ode in active_odes.items():
            out.append(ode(**dict_state))
        return out
    return simple_step_ode_function
