program Int2D;

uses
  Forms,
  xsect in 'xsect.pas' {Form1};

{$R *.RES}

begin
  Application.Initialize;
  Application.CreateForm(TForm1, Form1);
  Application.Run;
end.
