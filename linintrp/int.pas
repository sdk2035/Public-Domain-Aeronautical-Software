unit int;

(*   *)
(*   *)

interface

uses
  Windows, Messages, SysUtils, Classes, Graphics, Controls, Forms, Dialogs,
  StdCtrls, ExtCtrls;

type
  TForm1 = class(TForm)
    GroupBox1: TGroupBox;
    GroupBox2: TGroupBox;
    rgInterp: TRadioGroup;
    btnCompute: TButton;
    btnQuit: TButton;
    ebX1: TEdit;
    ebY1: TEdit;
    ebZ1: TEdit;
    ebX2: TEdit;
    ebY2: TEdit;
    ebZ2: TEdit;
    ebX3: TEdit;
    ebY3: TEdit;
    ebZ3: TEdit;
    Label1: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Label5: TLabel;
    Label6: TLabel;
    btnHelp: TButton;
    ebPercent: TEdit;
    procedure btnQuitClick(Sender: TObject);
    procedure btnComputeClick(Sender: TObject);
    procedure ebX3Change(Sender: TObject);
    procedure ebY3Change(Sender: TObject);
    procedure ebZ3Change(Sender: TObject);
    procedure ebPercentChange(Sender: TObject);
    procedure btnHelpClick(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  Form1: TForm1;

VAR
  x1,y1,z1 : EXTENDED;
  x2,y2,z2 : EXTENDED;
  x3,y3,z3 : EXTENDED;

implementation

{$R *.DFM}

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

procedure TForm1.btnQuitClick(Sender: TObject);
begin
  Application.Terminate
end;

procedure TForm1.btnComputeClick(Sender: TObject);
  VAR
    fraction,percent : EXTENDED;
begin
  DecodeEditBox(ebX1, x1, 'x1');
  DecodeEditBox(ebY1, y1, 'y1');
  DecodeEditBox(ebZ1, z1, 'z1');

  DecodeEditBox(ebX2, x2, 'x2');
  DecodeEditBox(ebY2, y2, 'y2');
  DecodeEditBox(ebZ2, z2, 'z2');

  CASE rgInterp.ItemIndex OF
    0 : BEGIN
         DecodeEditBox(ebX3, x3, 'x3');
         fraction:=(x3-x1)/(x2-x1);
         percent:=100.0*fraction;
         y3:=y1+fraction*(y2-y1);
         ebY3.Text:=FloatToStrF(y3,ffFixed,12,4);
         z3:=z1+fraction*(z2-z1);
         ebZ3.Text:=FloatToStrF(z3,ffFixed,12,4);
         ebPercent.Text:=FloatToStrF(percent,ffFixed,12,4);
         rgInterp.ItemIndex:=0;
         ebX3.SetFocus
        END;
    1 : BEGIN
         DecodeEditBox(ebY3, y3, 'y3');
         fraction:=(y3-y1)/(y2-y1);
         percent:=100.0*fraction;
         x3:=x1+fraction*(x2-x1);
         ebX3.Text:=FloatToStrF(x3,ffFixed,12,4);
         z3:=z1+fraction*(z2-z1);
         ebZ3.Text:=FloatToStrF(z3,ffFixed,12,4);
         ebPercent.Text:=FloatToStrF(percent,ffFixed,12,4);
         rgInterp.ItemIndex:=1;
         ebY3.SetFocus
        END;
    2 : BEGIN
         DecodeEditBox(ebZ3, z3, 'z3');
         fraction:=(z3-z1)/(z2-z1);
         percent:=100.0*fraction;
         x3:=x1+fraction*(x2-x1);
         ebX3.Text:=FloatToStrF(x3,ffFixed,12,4);
         y3:=y1+fraction*(y2-y1);
         ebY3.Text:=FloatToStrF(y3,ffFixed,12,4);
         ebPercent.Text:=FloatToStrF(percent,ffFixed,12,4);
         rgInterp.ItemIndex:=2;
         ebZ3.SetFocus
        END;
    3 : BEGIN
         DecodeEditBox(ebPercent,percent,'%');
         fraction:=0.01*percent;
         x3:=x1+fraction*(x2-x1);
         ebX3.Text:=FloatToStrF(x3,ffFixed,12,4);
         y3:=y1+fraction*(y2-y1);
         ebY3.Text:=FloatToStrF(y3,ffFixed,12,4);
         z3:=z1+fraction*(z2-z1);
         ebZ3.Text:=FloatToStrF(z3,ffFixed,12,4);
         rgInterp.ItemIndex:=3;
         ebPercent.SetFocus
        END;
  END;

end;

procedure TForm1.ebX3Change(Sender: TObject);
begin
  rgInterp.ItemIndex:=0
end;

procedure TForm1.ebY3Change(Sender: TObject);
begin
  rgInterp.ItemIndex:=1
end;

procedure TForm1.ebZ3Change(Sender: TObject);
begin
  rgInterp.ItemIndex:=2
end;

procedure TForm1.ebPercentChange(Sender: TObject);
begin
  rgInterp.ItemIndex:=3
end;

procedure TForm1.btnHelpClick(Sender: TObject);
begin
  MessageDlg('This tool helps you find points on a given line segment ' +
   'by interpolation. The line in 3D is defined by (x,y,z) of two points. ' +
   'On the left of the form are six boxes for the input of these numbers. ',
   mtInformation,[mbOK],0)
end;

end.
