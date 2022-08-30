import numpy as np
import plotly.io as pio
import plotly.graph_objects as go
import matplotlib as mpl
from dash import Dash, html, dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output
import dash_dangerously_set_inner_html


### You can find this application hosted on a web server
### https://population-app-dashboard.herokuapp.com/
### running this script locally allows the script to carry out computationally intensive tasks


########## styling
def colorFader(c1,c2,mix=0.0): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1=np.array(mpl.colors.to_rgb(c1))
    c2=np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)
################################

#Wright Fischer sampling
def extract(p, size):
    N=size
    randomlist = np.random.rand(2*N)
    newgen_alleles = [1 if x <= p else 0 for x in randomlist]
    # calcolo dei genotipi
    newgentype = [(newgen_alleles[i]+newgen_alleles[i+1]) for i in range(0, (2*N), 2)]
    #calcolo della nuova frequenza
    newp = sum(newgen_alleles)/(2*N)
    return newgen_alleles, newgentype, newp

#Frequency calculation at each generation
#each sample is a vector
def calcp_vector(p, gens, size,p0):
    pv = []
    pv.append(p0)
    for g in range(gens):
        a, b, p = extract(p, size)
        pv.append(p)
        if p == 0 or p == 1: break
    while len(pv) < gens+1:
        pv.append(pv[-1])
    return pv

#merge all the samples in a matrix
def calcp_matrix(size, p0, rep, gens):
    pm = []
    for r in range(rep+1):
        pv = calcp_vector(p0, gens,size,p0)
        pm.append(pv)
    return pm


#frequency vs generations plot
def calc_xx(gens):
    xx = np.arange(gens+1)
    return xx

def graph_freq(pm,size,rep, gens):
    xx = calc_xx(gens)
    fig = go.Figure()
    fig.update_layout(
        yaxis_range=[0, 1],
        template='plotly_white+grad',
        title=dict(
            text='Population Size N=%s' %size,
            font=dict(size=18, color='#5A5A5A'),
        )
    )
    fig.update_xaxes(
        title_text='Generations'
    )
    fig.update_yaxes(
        title_text='Allele Frequency (p<sub>A</sub>)'
    )

    for i in range(rep):
        fig.add_trace(go.Scatter(
            x=xx,
            y=pm[i],
            showlegend=False
        ))
    return fig

#Selection : deltaP calculation
def deltaP(p_vec, w_vec, w11, w12, w22):
    qq = 1 - p_vec[-1]
    w_mean = (p_vec[-1]) ** 2 * w11 + 2 * p_vec[-1] * qq * w12 + qq ** 2 * w22
    w_vec.append(w_mean)
    deltap = (p_vec[-1] * qq * (
                (p_vec[-1] * w11 + qq * w12) - (p_vec[-1] * w12 + qq * w22))) / w_mean
    return deltap

#Selection : frequency p calculation at each generation
def calcp_vector_sel(size, p, gens, w11, w12, w22):
    pv = []
    pv.append(p)
    w_vec = []
    for g in range(gens):
        a, b, p = extract(p, size)
        if p == 0 or p == 1:
            pv.append(p)
            break
        deltap = deltaP(pv,w_vec, w11, w12, w22)
        p = p + deltap
        pv.append(p)
    while len(pv) < gens+1:
        pv.append(pv[-1])
    return pv


def calcp_matrix_sel(size, w11, w12, w22, p0, rep, gens):
    pm = []
    for r in range(rep+1):
        pv = calcp_vector_sel(size, p0, gens, w11, w12, w22)
        pm.append(pv)
    return pm




#############################################################################
###########################################################################
#### The following part only contains the HTML structure to create the dashboard
#### All the functions needed for the genetic drift simulator are defined above


### Dashboard
dbc_css = "https://cdn.jsdelivr.net/gh/AnnMarieW/dash-bootstrap-templates/dbc.min.css"

app = Dash(__name__, external_stylesheets=[dbc.themes.MINTY, dbc_css])
server = app.server


app.layout = html.Div(
    [
        dbc.Row(
            [
                dbc.Alert(
                    [
                        html.H1('''Genetic Drift Simulation'''),
                        html.P(dash_dangerously_set_inner_html.DangerouslySetInnerHTML('''
                            This is a basic web app to explore the effects of genetic drifts on populations.
                            </br>The purpose of this simulation is entirely educative and does not intend to be an exhaustive description of evolution forces or any aspects of 
                            evolution in general. It is rather a user-friendly simple tool, useful to try out different scenarios. 
                            </br>The user can set all the variables needed to start the simulation.
                            You can also decide to integrate the effects of natural selections along with genetic drift.</br>

                            ''')
                        ),
                    ],
                    color='light',
                    style={'padding-left':'5%', 'padding-right':'5%', 'font-size':'14pt'},
                ),
            ],
            id='row_wright',
            style={'margin':'3%'}
        ),
        dbc.Row(
            [
                dbc.Col(
                    [html.A('Go to Plot', 
                        href='#first_div',
                        style={
                            'color':'#ffffff',
                            'padding':'2%',
                            'text-decoration':'none',
                            
                            'background-color':'#B0D5CB',
                            'border-radius':'10px',

                            }),

                    ],
                    width=6
                ),
                dbc.Col(html.Div([
                    dbc.Alert('This simulation is not intended for large-scale usage, thus it can only make use of a very limited computational power. Maxing out all the parameters may force the simulation to shut down.', color='primary')
                ]),
                    width=6
                )
            ],
            id='row_disclaimer',
            style={'margin':'1%'}
        ),
        dbc.Row(
            [
                dbc.Col(
                    [html.Div([
                    html.H2('Population Size'),
                    html.P(dash_dangerously_set_inner_html.DangerouslySetInnerHTML('''
                    Population size affects the way allele frequency varies between generations. </br>
                    Small populations are more affected by genetic drift than larger ones. This means that alleles are more likely to be lost or to achieve fixations (that is the case when the allele frequency reaches 1).</br>
                    If genetic drift is the only evolutionary force acting on the population, the fixation probability is equal to the starting frequency of the allele in the population.</br>
                    Since large populations are not really influenced by genetic drift, it is expected that in such situations the mean frequency for all the replicas throughout the generations is constant and approximately
                    equal to the starting frequency for the focus allele.
                    '''),
                    )
                ])],
                    width=6
                ),
                dbc.Col(
                    [html.Div(
                        [
                            dbc.Alert(html.Div(['Drag the slider or click on numbers to set population size. Note that it is log-scaled.'], style={'font-size':'14pt'}),color='light'),
                            dcc.Slider(0, 4, 0.01,
                            value=1,
                            marks={i:'{}'.format(10**i) for i in range(5)},
                            className='dbc',
                            id='Size'),
                        ],
                        style={'width':'80%'}
                        ), 
                        ],
                    width=6
                )
            ],
            id='row_size',
            style={'margin':'1%'}
        ),
        html.Hr(
            style={
                'background':'#91B3A1',
                'margin':'1%',
                'width':'80%',
                'height':'5px',
                'position':'relative',
                'margin':'auto'
                }
        ),
        dbc.Row(
            [
                dbc.Col(
                    [html.Div([
                    html.H2('Replicas Number'),
                    html.P(''' You can set the number of independent simulations you want to run simultaneously. These will be plotted together at the bottom of this page.
                    Simulations results are stochastic, this means you might need to replicate the simulation many times to observe the expected result (on average). '''),
                    html.P('High replicas number will slow the simulation.')
                ])],
                    width=6
                ),
                dbc.Col(
                    [html.Div(
                        [
                            dbc.Alert(html.Div(['Drag the slider to set the replicas number'], style={'font-size':'14pt'}),color='light'),
                            dcc.Slider(1,100,1,
                            value=5,
                            marks=None,
                            tooltip={'placement':'bottom','always_visible':True},
                            className='dbc',
                            id='Replicas'),
                        ],
                        style={'width':'80%'}
                        ), 
                        ],
                    width=6
                )
            ],
            id='row_replicas',
            style={'margin':'1%'}
        ),
        html.Hr(
            style={
                'background':'#91B3A1',
                'margin':'1%',
                'width':'80%',
                'height':'5px',
                'position':'relative',
                'margin':'auto'
                }
        ),
        dbc.Row(
            [
                dbc.Col(
                    [html.Div([
                    html.H2('Generations'),
                    html.P(dash_dangerously_set_inner_html.DangerouslySetInnerHTML(''' 
                        You can choose how long the simulation will run for. The time variable is presented in form of generations. 
                        </br>It is assumed that generations are non-overlapping and that all specimens can mate to give offspring.</br>
                    ''')
                    
                    ),
                ])],
                    width=6
                ),
                dbc.Col(
                    [html.Div(
                        [
                            dbc.Alert(html.Div(['Drag the slider or click on numbers to set the span of generations.'], style={'font-size':'14pt'}),color='light'),
                            dcc.Slider(1,300,1,
                            value=50,
                            marks={10:'10',50:'50',100:'100', **{i:'{}'.format(i) for i in range(200,400, 100)}},
                            tooltip={'placement':'bottom','always_visible':False},
                            id='Generations',
                            className='dbc'),
                        ],
                        style={'width':'80%'}
                        ), 
                        ],
                    width=6
                )
            ],
            id='row_generations',
            style={'margin':'1%'}
        ),
        html.Hr(
            style={
                'background':'#91B3A1',
                'margin':'1%',
                'width':'80%',
                'height':'5px',
                'position':'relative',
                'margin':'auto'
                }
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.Div([
                            html.H2('Drift and Selection'),
                            html.P(dash_dangerously_set_inner_html.DangerouslySetInnerHTML('''
                                In this simulation you can analyze genetic drift alone or coupled with natural selection.<br/>
                                If the focus allele has no adaptive behavior then selection should not influence its evolution. However, if you want, you can set different 
                                values of relative fitness for the three possible genotypes resulting from a biallelic system (AA, AB and BB).</br>
                                In this way you can observe how selection can be altered by genetic drift.
                            ''')
                            )
                        ]),
                    ],
                    width=6,
                ),
                dbc.Col(
                    [
                    dcc.Dropdown(
                        id='String',
                        options=[{
                            'label':'Genetic Drift only',
                            'value':'drift'
                        },
                        {
                            'label':'Genetic Drift with Selection',
                            'value':'selection'
                        }],
                        value='drift',
                        )
                    ],
                    style={'display':'inline-block'}
                    ),
                ],
                id='row_drift',
                style={'margin':'1%'}
        ),
        html.Hr(
            style={
                'background':'#91B3A1',
                'margin':'1%',
                'width':'80%',
                'height':'5px',
                'position':'relative',
                'margin':'auto'
                }
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.Div([
                            html.H2('Initial Allele Frequency'),
                            dash_dangerously_set_inner_html.DangerouslySetInnerHTML(['''You can change the starting frequency of focus allele A (p<sub>A</sub>)
                            </br>Assuming there are only 2 possible alleles for each locus the frequency for allele B is obtained by calculating 1-p<sub>A</sub>'''
                            ]),
                        ]),
                    ],
                    width=6,
                ),
                dbc.Col(
                    [
                    html.Div(
                    [
                    html.Div(['Allele A Frequency (p',html.Sub('A'),')']),
                    dcc.Slider(
                        0,1,0.01,
                        id='Frequency',
                        value=0.5,
                        tooltip={'placement':'bottom','always_visible':False},
                        marks={i/10 :'{}'.format(i/10) for i in range(0,11)},
                        className='dbc'
                        ),
                    ]
                    )
                    ],
                    ),
                ],
                id='row_init_freq',
                style={'margin':'1%'}
        ),
        html.Hr(
            style={
                'background':'#91B3A1',
                'margin':'1%',
                'width':'80%',
                'height':'5px',
                'position':'relative',
                'margin':'auto'
                }
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.Div([
                            html.H2('Genotype Relative Fitness'),
                            html.P([dash_dangerously_set_inner_html.DangerouslySetInnerHTML('''
                                The absolute fitness W of a genotype is defined as the proportional change in the abundance of that genotype over one generation attributable to selection. 
                                It can be thought as the average rate of growth per-capita.</br>
                                The Relative fitness w of a genotype in this case is defined as the ratio between its absolute fitness and the highest absolute fitness among the genotypes in the population.</br>
                                Thus at least one genotype should have relative fitness equal to 1. However, you can select any value for the relative fitness of the 3 genotypes of this simplified scenario
                                and they will be internally rescaled in order to respect the previous constraint [meaning that for the genotype triplet (AA,AB,BB) the same results will be observed (on average) using respectively the relative fitnesses of (1, 0.6, 0.4) or (0.5, 0.3, 0.2) ].
                            ''')]),
                        ]),
                    ],
                    width=6,
                ),
                dbc.Col(
                    dbc.Row(
                    [
                    html.Div(
                    [
                    html.Div('Genotype AA'),
                    dcc.Slider(
                        0,1,0.01,
                        tooltip={'placement':'bottom','always_visible':False},
                        id='FitnessAA',
                        marks=None,
                        value=0.3,
                        className='dbc'
                    )],
                    style={'width':'60%', 'margin':'1%',}),
                    html.Div([
                    html.Div('Genotype AB'),
                    dcc.Slider(
                        0,1,0.01,
                        tooltip={'placement':'bottom','always_visible':False},
                        id='FitnessAB',
                        marks=None,
                        value=0.3,
                        className='dbc'
                    )],
                    style={'width':'60%', 'margin':'1%', }),
                    html.Div([
                    html.Div('Genotype BB'),
                    dcc.Slider(
                        0,1,0.01,
                        tooltip={'placement':'bottom','always_visible':False},
                        id='FitnessBB',
                        marks=None,
                        value=0.3,
                        className='dbc'
                    )],
                    style={'width':'60%', 'margin':'1%',}),
                    ],
                    ),),
                ],
                id='row_fitness',
                style={'margin':'1%'}
        ),
        html.Div(id='first_div'),
    ],
    style={
            'margin':'auto',
            'width':'90vw',
            'padding-top':'3vh',
            'scroll-behavior':'smooth',
        }


)

def transform_log_value(v):
    return 10 ** v

@app.callback(
    Output('first_div','children'),
    Input('Size','value'),
    Input('String','value'),
    Input('Frequency','value'),
    Input('FitnessAA','value'),
    Input('FitnessAB','value'),
    Input('FitnessBB','value'),
    Input('Replicas','value'),
    Input('Generations','value')
)
def plot_drift(size, string, p0, w11, w12, w22, rep, gens):
    size=int(transform_log_value(size))
    if None in [size, string, p0, w11, w12, w22, rep]:
        return ''
    if size > 10000:
        return dbc.Alert('Please choose a smaller population size (Max 10000)', color='primary')

    set_pio_colorscale(rep)

    if string == 'drift':
        pm = calcp_matrix(size, p0,rep, gens)
    if string == 'selection':
        pm = calcp_matrix_sel(size, w11, w12, w22, p0, rep, gens)

    return dcc.Graph(figure=graph_freq(pm, size, rep, gens))
    
def set_pio_colorscale(rep):
    c1 = '#5a7976'
    c2 = '#b0dadd'
    pio.templates['grad'] = go.layout.Template(
        layout=dict(
            colorway=[colorFader(c1, c2, x/rep) for x in range(rep)]
        )
    )

@app.callback(
    Output('row_fitness','style'),
    Input('String','value')
)
def show_fitness(string):
    if string == 'selection':
        return {'visibility':'visible', 'margin':'1%'}
    else:
        return {'visibility':'hidden','height':'0.5vh'}


if __name__ == '__main__':
    app.run_server(debug=True)