object Form1: TForm1
  Left = 200
  Top = 100
  Width = 419
  Height = 363
  Caption = 'Intersection of two lines (2-D)'
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -13
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 120
  TextHeight = 16
  object TLabel
    Left = 216
    Top = 280
    Width = 3
    Height = 16
  end
  object TLabel
    Left = 216
    Top = 304
    Width = 3
    Height = 16
  end
  object Label9: TLabel
    Left = 275
    Top = 288
    Width = 3
    Height = 16
    Caption = ' '
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'MS Sans Serif'
    Font.Style = []
    ParentFont = False
  end
  object Label10: TLabel
    Left = 275
    Top = 312
    Width = 3
    Height = 16
    Caption = ' '
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'MS Sans Serif'
    Font.Style = []
    ParentFont = False
  end
  object GroupBox1: TGroupBox
    Left = 16
    Top = 8
    Width = 377
    Height = 129
    Caption = 'Line 1'
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'MS Sans Serif'
    Font.Style = []
    ParentFont = False
    TabOrder = 0
    object GroupBox3: TGroupBox
      Left = 24
      Top = 24
      Width = 153
      Height = 89
      Caption = 'First Point'
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'MS Sans Serif'
      Font.Style = []
      ParentFont = False
      TabOrder = 0
      object Label1: TLabel
        Left = 8
        Top = 24
        Width = 6
        Height = 16
        Caption = 'x'
      end
      object Label2: TLabel
        Left = 8
        Top = 60
        Width = 7
        Height = 16
        Caption = 'y'
      end
      object ebL1P1x: TEdit
        Left = 24
        Top = 24
        Width = 121
        Height = 24
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'MS Sans Serif'
        Font.Style = []
        ParentFont = False
        TabOrder = 0
        OnChange = ebL1P1xChange
      end
      object ebL1P1y: TEdit
        Left = 24
        Top = 56
        Width = 121
        Height = 24
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'MS Sans Serif'
        Font.Style = []
        ParentFont = False
        TabOrder = 1
        OnChange = ebL1P1yChange
      end
    end
    object GroupBox4: TGroupBox
      Left = 200
      Top = 24
      Width = 153
      Height = 89
      Caption = 'Second Point'
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'MS Sans Serif'
      Font.Style = []
      ParentFont = False
      TabOrder = 1
      object Label3: TLabel
        Left = 8
        Top = 24
        Width = 6
        Height = 16
        Caption = 'x'
      end
      object Label4: TLabel
        Left = 8
        Top = 60
        Width = 7
        Height = 16
        Caption = 'y'
      end
      object ebL1P2x: TEdit
        Left = 24
        Top = 24
        Width = 121
        Height = 24
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'MS Sans Serif'
        Font.Style = []
        ParentFont = False
        TabOrder = 0
        OnChange = ebL1P2xChange
      end
      object ebL1P2y: TEdit
        Left = 24
        Top = 56
        Width = 121
        Height = 24
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'MS Sans Serif'
        Font.Style = []
        ParentFont = False
        TabOrder = 1
        OnChange = ebL1P2yChange
      end
    end
  end
  object GroupBox2: TGroupBox
    Left = 16
    Top = 152
    Width = 377
    Height = 129
    Caption = 'Line 2'
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'MS Sans Serif'
    Font.Style = []
    ParentFont = False
    TabOrder = 1
    object GroupBox5: TGroupBox
      Left = 24
      Top = 16
      Width = 153
      Height = 89
      Caption = 'First Point'
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'MS Sans Serif'
      Font.Style = []
      ParentFont = False
      TabOrder = 0
      object Label5: TLabel
        Left = 8
        Top = 24
        Width = 6
        Height = 16
        Caption = 'x'
      end
      object Label6: TLabel
        Left = 8
        Top = 60
        Width = 7
        Height = 16
        Caption = 'y'
      end
      object ebL2P1x: TEdit
        Left = 24
        Top = 24
        Width = 121
        Height = 24
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'MS Sans Serif'
        Font.Style = []
        ParentFont = False
        TabOrder = 0
        OnChange = ebL2P1xChange
      end
      object ebL2P1y: TEdit
        Left = 24
        Top = 56
        Width = 121
        Height = 24
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'MS Sans Serif'
        Font.Style = []
        ParentFont = False
        TabOrder = 1
        OnChange = ebL2P1yChange
      end
    end
    object GroupBox6: TGroupBox
      Left = 200
      Top = 16
      Width = 153
      Height = 89
      Caption = 'Second Point'
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'MS Sans Serif'
      Font.Style = []
      ParentFont = False
      TabOrder = 1
      object Label7: TLabel
        Left = 8
        Top = 24
        Width = 6
        Height = 16
        Caption = 'x'
      end
      object Label8: TLabel
        Left = 8
        Top = 60
        Width = 7
        Height = 16
        Caption = 'y'
      end
      object ebL2P2x: TEdit
        Left = 24
        Top = 24
        Width = 121
        Height = 24
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'MS Sans Serif'
        Font.Style = []
        ParentFont = False
        TabOrder = 0
        OnChange = ebL2P2xChange
      end
      object ebL2P2y: TEdit
        Left = 24
        Top = 56
        Width = 121
        Height = 24
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'MS Sans Serif'
        Font.Style = []
        ParentFont = False
        TabOrder = 1
        OnChange = ebL2P2yChange
      end
    end
  end
  object btnCompute: TButton
    Left = 16
    Top = 296
    Width = 75
    Height = 25
    Hint = 'Compute the intersection'
    Caption = 'Compute'
    Default = True
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'MS Sans Serif'
    Font.Style = []
    ParentFont = False
    ParentShowHint = False
    ShowHint = True
    TabOrder = 2
    OnClick = btnComputeClick
  end
  object btnQuit: TButton
    Left = 104
    Top = 296
    Width = 65
    Height = 25
    Hint = 'Quit the program'
    Cancel = True
    Caption = 'Quit'
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'MS Sans Serif'
    Font.Style = []
    ParentFont = False
    ParentShowHint = False
    ShowHint = True
    TabOrder = 3
    OnClick = btnQuitClick
  end
  object btnHelp: TButton
    Left = 184
    Top = 296
    Width = 65
    Height = 25
    Hint = 'Show instructions on the use of the program'
    Caption = 'Help'
    ParentShowHint = False
    ShowHint = True
    TabOrder = 4
    OnClick = btnHelpClick
  end
end
