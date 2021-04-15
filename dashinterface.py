# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_table
import plotly.express as px
from dash.dependencies import Input, Output

import os
from flask import Flask


import urllib.parse

import json
import pandas as pd

from molmass import Formula
import IsoSpecPy

from Molecule import Molecule, molecular_factory_dict
from adducts import get_adduct_mass

from app import app

dash_app = dash.Dash(__name__, server=app, external_stylesheets=[dbc.themes.BOOTSTRAP], url_base_pathname='/dashinterface/')
dash_app.title = 'Calculator'

NAVBAR = dbc.Navbar(
    children=[
        dbc.NavbarBrand(
            html.Img(src="https://gnps-cytoscape.ucsd.edu/static/img/GNPS_logo.png", width="120px"),
            href="https://gnps.ucsd.edu"
        ),
        dbc.Nav(
            [
                dbc.NavItem(dbc.NavLink("GNPS Mass Spec Calculator", href="#")),
            ],
        navbar=True)
    ],
    color="light",
    dark=False,
    sticky="top",
)

DASHBOARD = [
    dbc.CardHeader(html.H5("GNPS Mass Spec Calculator")),
    dbc.CardBody(
        [
            html.Div(id='version', children="Version - Release_2"),
            html.Br(),
            dbc.InputGroup(
                [
                    dbc.InputGroupAddon("Molecular Formula", addon_type="prepend"),
                    dbc.Input(placeholder="Enter Formula", id="formula_entry"),
                ],
                className="mb-3",
            ),
            html.Br(),
            dbc.InputGroup(
                [
                    dbc.InputGroupAddon("SMILES", addon_type="prepend"),
                    dbc.Input(placeholder="Enter SMILES structure", id="smiles_entry", value=""),
                ],
                className="mb-3",
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupAddon("InChI", addon_type="prepend"),
                    dbc.Input(placeholder="Enter InChI structure", id="inchi_entry", value=""),
                ],
                className="mb-3",
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupAddon("InChIKey", addon_type="prepend"),
                    dbc.Input(placeholder="Enter InChIKey structure", id="inchikey_entry", value=""),
                ],
                className="mb-3",
            ),
            
            html.Hr(),

            html.H4("Stucture Information"),
            html.Hr(),
            dcc.Loading(
                className="mb-3",
                id="structureinfo",
                children=[html.Div([html.Div(id="loading-output-44")])],
                type="default",
            ),

            html.Hr(),
            html.H4("Monoisotopic Adduct Table"),
            html.Hr(),
            dcc.Loading(
                className="mb-3",
                id="massspecinfo",
                children=[html.Div([html.Div(id="loading-output-3")])],
                type="default",
            ),

            html.Hr(),
            html.H4("Isotopologue Table"),
            html.Hr(),
            dcc.Slider(
                min=1,
                max=6,
                value=3,
                step=0.1,
                id="resolution-slider",
                marks={
                    1: {'label': '10'},
                    2: {'label': '100'},
                    3: {'label': '1000'},
                    4: {'label': '10000'},
                    5: {'label': '100000'},
                    6: {'label': '1000000'},
                },

            ),
            html.Br(),
            dcc.Loading(
                className="mb-3",
                id="isotopologueinfo",
                children=[html.Div([html.Div(id="loading-output-4")])],
                type="default",
            )
        ]
    )
]

BODY = dbc.Container(
    [
        dbc.Row([dbc.Col(dbc.Card(DASHBOARD)),], style={"marginTop": 30}),
    ],
    className="mt-12",
)

dash_app.layout = html.Div(children=[NAVBAR, BODY])


@dash_app.callback(
    [Output('structureinfo', 'children')],
    [Input('smiles_entry', 'value'), Input('inchi_entry', 'value'), Input('inchikey_entry', 'value')]
)
def generate_structure_information(smiles_entry, inchi_entry, inchikey_entry):
    structure_dict = {}
    if len(smiles_entry) > 0:
        structure_dict["smiles"] = smiles_entry
    if len(inchi_entry) > 0:
        structure_dict["inchi"] = inchi_entry
    if len(inchikey_entry) > 0:
        structure_dict["inchikey"] = inchikey_entry

    m = molecular_factory_dict(structure_dict)

    result_list = []
    result_list.append({"structure" : "InChI", "value" : m.inchi})
    result_list.append({"structure" : "SMILES", "value" : m.smiles})
    result_list.append({"structure" : "InChIKey", "value" : m.inchikey})

    structure_df = pd.DataFrame(result_list)

    table = dbc.Table.from_dataframe(structure_df, striped=True, bordered=True, hover=True)
    image_src = "/structureimg?{}".format(urllib.parse.urlencode(structure_dict))
    if len(structure_dict) > 0:
        img_html = html.Img(src=image_src)
    else:
        img_html = "Please Enter Structure"

    return [
        dbc.Row([
            dbc.Col(table),
            dbc.Col(img_html)
        ])
        
    ]



# This function will generate monoisotopic masses for a set of adducts
@dash_app.callback(
    [Output('massspecinfo', 'children')],
    [Input('formula_entry', 'value'), Input('smiles_entry', 'value'), Input('inchi_entry', 'value'), Input('inchikey_entry', 'value')],
)
def generate_adduct_information(formula_entry, smiles_entry, inchi_entry, inchikey_entry):
    exact_mass = 0

    if formula_entry is not None and len(formula_entry):
        f = Formula(formula_entry)
        exact_mass = f.isotope.mass
    else:
        # Getting exact mass
        structure_dict = {}
        if len(smiles_entry) > 0:
            structure_dict["smiles"] = smiles_entry
        if len(inchi_entry) > 0:
            structure_dict["inchi"] = inchi_entry
        if len(inchikey_entry) > 0:
            structure_dict["inchikey"] = inchikey_entry
        m = molecular_factory_dict(structure_dict)
        exact_mass = float(m.exact_mass)

    adducts_to_report = ["M", "M+H", "M+Na", "M+K", "M+NH4", "M-H", "M+Br", "M+Cl"]
    output_list = []

    for adduct in adducts_to_report:
        adduct_mass, charge = get_adduct_mass(exact_mass, adduct)
        output_dict = {}
        output_dict["adduct"] = adduct
        output_dict["charge"] = charge
        output_dict["mz"] = adduct_mass
        output_list.append(output_dict)

    table_fig = dash_table.DataTable(
        columns=[
            {"name": i, "id": i, "deletable": True, "selectable": True} for i in ["adduct", "charge", "mz"]
        ],
        data=output_list,
        editable=True,
        filter_action="native",
        sort_action="native",
        sort_mode="multi",
        column_selectable="single",
        selected_columns=[],
        selected_rows=[],
        page_action="native",
        page_current= 0,
        page_size= 10,
    )

    return [table_fig]


@dash_app.callback(
    [Output('isotopologueinfo', 'children')],
    [
        Input('formula_entry', 'value'), 
        Input('smiles_entry', 'value'), 
        Input('inchi_entry', 'value'), 
        Input('inchikey_entry', 'value'), 
        Input("resolution-slider", "value")
    ],
)
def generate_isotopologues(formula_entry, smiles_entry, inchi_entry, inchikey_entry, resolution_entry):
    formula = ""
    if formula_entry is not None and len(formula_entry):
        formula = formula_entry
    else:
        # Getting exact mass
        structure_dict = {}
        if len(smiles_entry) > 0:
            structure_dict["smiles"] = smiles_entry
        if len(inchi_entry) > 0:
            structure_dict["inchi"] = inchi_entry
        if len(inchikey_entry) > 0:
            structure_dict["inchikey"] = inchikey_entry
        m = molecular_factory_dict(structure_dict)
        formula = (m.formula)

    i = IsoSpecPy.IsoTotalProb(formula = formula, # The formula for glucose, sans the radiolabel atoms                            # And the rest of parameters for configuration
                            prob_to_cover = 0.99, 
                            get_confs=True)
    output_list = []
    for mass, prob, conf in i:
        output_dict = {}
        output_dict["prob"] = prob
        output_dict["mz"] = mass - 0.00054858
        output_list.append(output_dict)

    table_fig = dash_table.DataTable(
        columns=[
            {"name": i, "id": i, "deletable": True, "selectable": True} for i in ["mz", "prob"]
        ],
        data=output_list,
        editable=True,
        filter_action="native",
        sort_action="native",
        sort_mode="multi",
        column_selectable="single",
        selected_columns=[],
        selected_rows=[],
        page_action="native",
        page_current= 0,
        page_size= 10,
    )

    # Drawing Figure
    main_mz = output_list[0]["mz"]
    true_resolution = (10 ** float(resolution_entry))
    delta_m = main_mz / true_resolution
    sigma = delta_m/2.355

    display_bins = 0.02
    display_bins = sigma

    import numpy as np
    mz_grid = np.arange(output_list[0]["mz"] - 1,
                        output_list[-1]["mz"] + 1, display_bins)
    intensity = np.zeros_like(mz_grid)

    for peak in output_list:
        # Add gaussian peak shape centered around each theoretical peak
        intensity += peak["prob"] * np.exp(-(mz_grid - peak["mz"]) ** 2 / (2 * sigma)
                ) / (np.sqrt(2 * np.pi) * sigma)

    # Normalize profile to 0-100
    intensity = (intensity / intensity.max()) * 100

    df = pd.DataFrame()
    df["mz"] = mz_grid
    df["intensity"] = intensity

    line_fig = px.line(df, x="mz", y="intensity", title='Isotopologue Distribution - {} - Resolution - {}'.format(formula, true_resolution))

    return [["Resolution Entry - {}".format(true_resolution), html.Hr(), table_fig, dcc.Graph(figure=line_fig)]]

if __name__ == "__main__":
    app.run_server(debug=True, port=5000, host="0.0.0.0")
