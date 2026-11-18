from streamlit.testing.v1 import AppTest
from pathlib import Path
from statespace.examples import EXEMPLOS

def test_full_interface_examples_and_capture():
    app=AppTest.from_file(str(Path(__file__).resolve().parents[1]/'app.py'),default_timeout=40).run()
    assert not app.exception
    for name in EXEMPLOS:
        app.selectbox(key='example').select(name).run()
        app.button[0].click().run()
        app.button[1].click().run()
        assert not app.exception, (name,list(app.exception))
        assert not app.error, (name,[x.value for x in app.error])
        assert len(app.latex)>10
    app.checkbox(key='capture').check().run()
    assert not app.exception
    assert any(x.label=='← Voltar à edição' for x in app.button)
    next(x for x in app.button if x.label=='← Voltar à edição').click().run()
    assert not app.checkbox(key='capture').value
