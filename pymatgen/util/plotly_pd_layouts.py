# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Default plotly layouts for Binary (2), Ternary (3), and Quaternary (4) phase
diagrams.
"""

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
        orientation="v",
        x=0.1,
        y=0.99,
        traceorder="reversed",
        xanchor="left",
        yanchor="top",
    ),
    scene_camera=dict(
        center=dict(x=-0.1, y=-0.1, z=0),
        eye=dict(x=-0.1, y=-0.1, z=1),
        projection=dict(type="orthographic"),
    ),
    scene=dict(xaxis=default_3d_axis, yaxis=default_3d_axis, zaxis=default_3d_axis),
    scene_aspectratio=dict(x=1.7, y=1.7, z=1.2),
)

colorscale = [
    [0.0, "#008d00"],
    [0.1111111111111111, "#4b9f3f"],
    [0.2222222222222222, "#73b255"],
    [0.3333333333333333, "#97c65b"],
    [0.4444444444444444, "#b9db53"],
    [0.5555555555555556, "#ffdcdf"],
    [0.6666666666666666, "#ffb8bf"],
    [0.7777777777777778, "#fd92a0"],
    [0.8888888888888888, "#f46b86"],
    [1.0, "#e24377"],
]

default_binary_marker_settings = dict(
    mode="markers",
    marker=dict(size=6, colorscale=colorscale, line=dict(width=4, color="black")),
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

empty_plot_style = {
    "xaxis": {"visible": False},
    "yaxis": {"visible": False},
    "paper_bgcolor": "rgba(0,0,0,0)",
    "plot_bgcolor": "rgba(0,0,0,0)",
}
