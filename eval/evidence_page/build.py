#!/usr/bin/env python3
"""Fill template.html with cases_trna.json (from gather.py + add_trna.py).  usage: build.py <cases_trna.json> <out.html>"""
import sys, json, os
d = json.load(open(sys.argv[1])); t = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'template.html')).read()
open(sys.argv[2], 'w').write(t.replace('/*DATA*/', json.dumps(d['cases'])).replace('/*TRNA*/', json.dumps(d['trna_summary'])))
