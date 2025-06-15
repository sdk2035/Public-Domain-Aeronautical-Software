object Form1: TForm1
  Left = 192
  Top = 107
  Width = 640
  Height = 480
  Caption = 'Line Interpolation Utility Tool'
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -14
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 120
  TextHeight = 16
  object GroupBox1: TGroupBox
    Left = 39
    Top = 30
    Width = 228
    Height = 139
    Hint = 'defines (x,y,z) of the first point'
    Caption = 'Endpoint 1'
    ParentShowHint = False
    ShowHint = True
    TabOrder = 0
    object Label1: TLabel
      Left = 10
      Top = 20
      Width = 6
      Height = 16
      Caption = 'x'
    end
    object Label2: TLabel
      Left = 10
      Top = 59
      Width = 7
      Height = 16
      Caption = 'y'
    end
    object Label3: TLabel
      Left = 10
      Top = 98
      Width = 6
      Height = 16
      Caption = 'z'
    end
    object ebX1: TEdit
      Left = 30
      Top = 20
      Width = 148
      Height = 24
      TabOrder = 0
    end
    object ebY1: TEdit
      Left = 30
      Top = 59
      Width = 148
      Height = 24
      TabOrder = 1
    end
    object ebZ1: TEdit
      Left = 30
      Top = 98
      Width = 148
      Height = 24
      TabOrder = 2
    end
  end
  object GroupBox2: TGroupBox
    Left = 39
    Top = 207
    Width = 228
    Height = 139
    Hint = 'defines the (x,y,z) coordinates of the second point'
    Caption = 'Endpoint 2'
    ParentShowHint = False
    ShowHint = True
    TabOrder = 1
    object Label4: TLabel
      Left = 10
      Top = 20
      Width = 6
      Height = 16
      Caption = 'x'
    end
    object Label5: TLabel
      Left = 10
      Top = 59
      Width = 7
      Height = 16
      Caption = 'y'
    end
    object Label6: TLabel
      Left = 10
      Top = 98
      Width = 6
      Height = 16
      Caption = 'z'
    end
    object ebX2: TEdit
      Left = 30
      Top = 20
      Width = 148
      Height = 24
      TabOrder = 0
    end
    object ebY2: TEdit
      Left = 30
      Top = 59
      Width = 148
      Height = 24
      TabOrder = 1
    end
    object ebZ2: TEdit
      Left = 30
      Top = 98
      Width = 148
      Height = 24
      TabOrder = 2
    end
  end
  object rgInterp: TRadioGroup
    Left = 374
    Top = 39
    Width = 228
    Height = 170
    Caption = 'Interpolation Point'
    Items.Strings = (
      'x'
      'y'
      'z'
      '%')
    TabOrder = 2
  end
  object btnCompute: TButton
    Left = 374
    Top = 273
    Width = 92
    Height = 30
    Hint = 'Compute the interpolation point'
    Caption = 'Compute'
    Default = True
    ParentShowHint = False
    ShowHint = True
    TabOrder = 7
    OnClick = btnComputeClick
  end
  object btnQuit: TButton
    Left = 510
    Top = 274
    Width = 92
    Height = 31
    Hint = 'Quit the Application'
    Cancel = True
    Caption = 'Quit'
    ParentShowHint = False
    ShowHint = True
    TabOrder = 8
    OnClick = btnQuitClick
  end
  object ebX3: TEdit
    Left = 433
    Top = 59
    Width = 149
    Height = 24
    TabOrder = 3
    OnChange = ebX3Change
  end
  object ebY3: TEdit
    Left = 433
    Top = 98
    Width = 149
    Height = 24
    TabOrder = 4
    OnChange = ebY3Change
  end
  object ebZ3: TEdit
    Left = 433
    Top = 138
    Width = 149
    Height = 24
    TabOrder = 5
    OnChange = ebZ3Change
  end
  object btnHelp: TButton
    Left = 374
    Top = 315
    Width = 92
    Height = 31
    Hint = 'Display instructions for using the program'
    Caption = 'Help'
    ParentShowHint = False
    ShowHint = True
    TabOrder = 9
    OnClick = btnHelpClick
  end
  object ebPercent: TEdit
    Left = 432
    Top = 176
    Width = 153
    Height = 24
    TabOrder = 6
    OnChange = ebPercentChange
  end
end
