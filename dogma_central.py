from graphviz import Digraph

fluxo = Digraph('dogma', format='png')

fluxo.attr(rankdir='LR', bgcolor='lightgrey', fontname='Helvetica')

fluxo.node('DNA', shape='box', style='filled', color='lightblue', fontname='Helvetica-Bold')
fluxo.node('RNA', shape='ellipse', style='filled', color='lightgreen', fontname='Helvetica-Bold')
fluxo.node('Proteína', shape='diamond', style='filled', color='lightpink', fontname='Helvetica-Bold')

fluxo.edge('DNA', 'RNA', label='Transcrição')
fluxo.edge('RNA', 'Proteína', label='Tradução')

fluxo.render('dogma_central', format='png')
