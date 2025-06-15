unit xsect;

interface

uses
  Windows, Messages, SysUtils, Classes, Graphics, Controls, Forms, Dialogs,
  StdCtrls;

type
  TForm1 = class(TForm)
    GroupBox1: TGroupBox;
    GroupBox2: TGroupBox;
    GroupBox3: TGroupBox;
    GroupBox4: TGroupBox;
    GroupBox5: TGroupBox;
    GroupBox6: TGroupBox;
    ebL1P1x: TEdit;
    ebL1P1y: TEdit;
    ebL1P2x: TEdit;
    ebL1P2y: TEdit;
    ebL2P1x: TEdit;
    ebL2P1y: TEdit;
    ebL2P2x: TEdit;
    ebL2P2y: TEdit;
    btnCompute: TButton;
    btnQuit: TButton;
    Label1: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Label5: TLabel;
    Label6: TLabel;
    Label7: TLabel;
    Label8: TLabel;
    Label9: TLabel;
    Label10: TLabel;
    btnHelp: TButton;
    procedure btnQuitClick(Sender: TObject);
    procedure btnComputeClick(Sender: TObject);
    procedure ebL1P1xChange(Sender: TObject);
    procedure ebL1P1yChange(Sender: TObject);
    procedure ebL1P2xChange(Sender: TObject);
    procedure ebL1P2yChange(Sender: TObject);
    procedure ebL2P1xChange(Sender: TObject);
    procedure ebL2P1yChange(Sender: TObject);
    procedure ebL2P2xChange(Sender: TObject);
    procedure ebL2P2yChange(Sender: TObject);
    procedure btnHelpClick(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  Form1: TForm1;

implementation

{$R *.DFM}

VAR
  L1P1x,L1P1y : EXTENDED;
  L1P2x,L1P2y : EXTENDED;
  L2P1x,L2P1y : EXTENDED;
  L2P2x,L2P2y : EXTENDED;


PROCEDURE ErrorMsg(const s : STRING);
BEGIN
  MessageDlg(s, mtError, [mbOk], 0)
END;   (* ------------------------------------- End of Procedure Error *)

PROCEDURE DecodeEditBox(CONST eb : TEdit;
                        VAR a : EXTENDED;
                        CONST msg :STRING);
  VAR
    errCode : INTEGER;
BEGIN
  Val(eb.Text,a,errCode);
  IF errCode <> 0 THEN
    BEGIN
      ErrorMsg('invalid Entry for '+msg);
      eb.SetFocus;
      Exit
    END;
END;   (* -------------- End of Procedure DecodeErrorBoxes *)


{ The following procedure was inspired by the routine xlines.c
   by  Mukesh Prasad in Graphics Gems II. However, it looks for an 
   intersection of the two infinite lines rather than the line segments }
PROCEDURE LinesIntersect(const x1,y1,x2,y2: EXTENDED; { first line}
                         const x3,y3,x4,y4: EXTENDED; { second line }
                         VAR code : INTEGER; { =0 OK; =1 lines parallel}
                         VAR x,y : EXTENDED); { intersection point }

  VAR
    a1, a2, b1, b2, c1, c2 : EXTENDED; { Coefficients of line eqns.}
    denom : EXTENDED;

BEGIN
  a1:= y2-y1;
  b1:= x1-x2;
  c1:= x2*y1 - x1*y2;  { a1*x + b1*y + c1 = 0 is line 1 }

  a2:= y4-y3;
  b2:= x3-x4;
  c2:= x4*y3 - x3*y4;  { a2*x + b2*y + c2 = 0 is line 2 }

  denom:= a1*b2 - a2*b1;
  IF denom = 0 THEN
    BEGIN
      code:=1;
      EXIT
    END;

  x:=(b1*c2 - b2*c1)/denom;
  y:=(a2*c1 - a1*c2)/denom;
  code:=0
END;   (* ------ End Procedure LinesIntersect *)


procedure TForm1.btnQuitClick(Sender: TObject);
begin
  Application.Terminate
end;

PROCEDURE TForm1.btnComputeClick(Sender: TObject);
  VAR
    code : INTEGER;
    x,y : EXTENDED;
BEGIN
  DecodeEditBox(ebL1P1x,L1P1x,'line 1, point 1, x');
  DecodeEditBox(ebL1P1y,L1P1y,'line 1, point 1, y');
  DecodeEditBox(ebL1P2x,L1P2x,'line 1, point 2, x');
  DecodeEditBox(ebL1P2y,L1P2y,'line 1, point 2, y');
  DecodeEditBox(ebL2P1x,L2P1x,'line 2, point 1, x');
  DecodeEditBox(ebL2P1y,L2P1y,'line 2, point 1, y');
  DecodeEditBox(ebL2P2x,L2P2x,'line 2, point 2, x');
  DecodeEditBox(ebL2P2y,L2P2y,'line 2, point 2, y');

  LinesIntersect(L1P1x,L1P1y, L1P2x,L1P2y,
                 L2P1x,L2P1y, L2P2x,L2P2y, code, x,y);
  IF code <> 0 THEN
    ErrorMsg('The lines are parallel')
  ELSE
    BEGIN
      Label9.Caption:='x='+FloatToStrF(x,ffFixed,12,5);
      Label10.Caption:='y='+FloatToStrF(y,ffFixed,12,5)
    END
END;


procedure TForm1.ebL1P1xChange(Sender: TObject);
begin
  label9.Caption:=' '; label10.Caption:=' '
end;

procedure TForm1.ebL1P1yChange(Sender: TObject);
begin
  label9.Caption:=' '; label10.Caption:=' '
end;

procedure TForm1.ebL1P2xChange(Sender: TObject);
begin
  label9.Caption:=' '; label10.Caption:=' '
end;

procedure TForm1.ebL1P2yChange(Sender: TObject);
begin
  label9.Caption:=' '; label10.Caption:=' '
end;

procedure TForm1.ebL2P1xChange(Sender: TObject);
begin
  label9.Caption:=' '; label10.Caption:=' '
end;

procedure TForm1.ebL2P1yChange(Sender: TObject);
begin
  label9.Caption:=' '; label10.Caption:=' '
end;

procedure TForm1.ebL2P2xChange(Sender: TObject);
begin
  label9.Caption:=' '; label10.Caption:=' '
end;

procedure TForm1.ebL2P2yChange(Sender: TObject);
begin
  label9.Caption:=' '; label10.Caption:=' '
end;

procedure TForm1.btnHelpClick(Sender: TObject);
begin
  MessageDlg('This tool finds the intersection of two lines in 2D. ' +
   'Each line is defined by two points. ' +
   'Enter the (x,y) coordinates of each of the four points and ' +
   'click the compute button (or press Enter) to display the intersection. ' +
   'Click Quit (or press Esc) to quit the program. ' +
   'The Tab key will move the focus. ' +
   'An error is indicated if any of the boxes are unfilled or contain ' +
   'invalid characters. An error is indicated if the lines are parallel.'
     ,mtInformation,[mbOK],0)
end;

END.
