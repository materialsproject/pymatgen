# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Default plotly layouts for Binary (2), Ternary (3), and Quaternary (4) phase
diagrams.
"""

import numpy as np

default_binary_layout = dict(
    xaxis={
        "title": "Fraction",
        "anchor": "y",
        "mirror": "ticks",
        "nticks": 8,
        "showgrid": True,
        "showline": True,
        "side": "bottom",
        "tickfont": {"size": 16.0},
        "ticks": "inside",
        "titlefont": {"color": "#000000", "size": 20.0},
        "type": "linear",
        "zeroline": False,
        "gridcolor": "rgba(0,0,0,0.1)",
    },
    yaxis={
        "title": "Formation energy (eV/atom)",
        "anchor": "x",
        "mirror": "ticks",
        "showgrid": True,
        "showline": True,
        "side": "left",
        "tickfont": {"size": 16.0},
        "ticks": "inside",
        "titlefont": {"color": "#000000", "size": 20.0},
        "type": "linear",
        "gridcolor": "rgba(0,0,0,0.1)",
    },
    hovermode="closest",
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
    autosize=True,
    showlegend=True,
    legend=dict(
        orientation="h",
        traceorder="reversed",
        x=0.1,
        y=1.05,
        xanchor="right",
        tracegroupgap=5,
    ),
    margin=dict(b=10, l=10, pad=20, t=20, r=10),
)

default_3d_axis = dict(
    title=None,
    visible=False,
    autorange=True,
    showgrid=False,
    zeroline=False,
    showline=False,
    ticks="",
    showaxeslabels=False,
    showticklabels=False,
    showspikes=False,
)

default_ternary_layout = dict(
    autosize=True,
    hovermode="closest",
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
    margin=dict(b=0, l=0, pad=0, t=0, r=0),
    showlegend=True,
    legend=dict(
        orientation="h",
        x=0.5,
        y=0.0,
        traceorder="reversed",
        xanchor="center",
        yanchor="top",
    ),
    scene_camera=dict(
        center=dict(x=0, y=0, z=0),
        eye=dict(x=0, y=0, z=1),
        projection=dict(type="orthographic"),
    ),
    scene=dict(xaxis=default_3d_axis, yaxis=default_3d_axis, zaxis=default_3d_axis),
    scene_aspectratio=dict(x=1.7, y=1.7, z=1.4),
)

stable_colors = ['#0c8c00', '#d8ffd4', '#ffffff']
stable_colorscale = list(zip(np.linspace(0, 1, len(stable_colors)), stable_colors))

unstable_colors = ['#fcf1a9', '#ff813d', '#ff0000']
unstable_colorscale = list(zip(np.linspace(0, 1, len(unstable_colors)),
                               unstable_colors))

default_binary_marker_settings = dict(
    mode="markers",
    marker=dict(size=6, line=dict(width=4, color="black")),
    hoverinfo="text",
    hoverlabel=dict(font=dict(size=14)),
    showlegend=True,
)

default_ternary_marker_settings = default_binary_marker_settings.copy()

default_quaternary_marker_settings = default_binary_marker_settings.copy()
default_quaternary_marker_settings.update(dict(line=dict(width=2, color="black")))

default_quaternary_layout = dict(
    autosize=True,
    height=450,
    hovermode="closest",
    margin=dict(b=30, l=30, pad=0, t=0, r=20),
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
    showlegend=True,
    legend=dict(
        orientation="v",
        x=0.1,
        y=0.99,
        traceorder="reversed",
        xanchor="left",
        yanchor="top",
    ),
    scene=dict(xaxis=default_3d_axis, yaxis=default_3d_axis, zaxis=default_3d_axis),
)

default_annotation_layout = {
    "align": "center",
    "opacity": 0.7,
    "showarrow": False,
    "xanchor": "right",
    "yanchor": "auto",
    "xshift": -10,
    "yshift": -10,
    "xref": "x",
    "yref": "y",
}

empty_plot_style = {
    "xaxis": {"visible": False},
    "yaxis": {"visible": False},
    "paper_bgcolor": "rgba(0,0,0,0)",
    "plot_bgcolor": "rgba(0,0,0,0)",
}
