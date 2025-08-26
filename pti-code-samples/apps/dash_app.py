# Minimal Dash app to browse QC/DE (skeleton)
import dash
from dash import html, dcc, Input, Output
import pandas as pd

app = dash.Dash(__name__)
df = pd.DataFrame({'gene':['A','B','C'],'logFC':[2.1,-1.2,0.7],'pval':[1e-4,0.03,0.2]})

app.layout = html.Div([
    html.H3('DE Browser'),
    dcc.Input(id='q', placeholder='Filter gene...', type='text'),
    html.Div(id='tbl')
])

@app.callback(Output('tbl','children'), Input('q','value'))
def render(q):
    sub = df if not q else df[df.gene.str.contains(q, case=False, na=False)]
    return html.Pre(sub.to_string(index=False))

if __name__ == '__main__':
    app.run_server(debug=True)
