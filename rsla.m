function []=rsla(varargin)
    clc;
    
    if nargin
        if strcmp(varargin{1},'init')
            rslastarter
        elseif strcmp(varargin{1},'analyze')
            rslaanalyzer()
        end
    else
        WelcomeText='Welcome to the resting-state lag-analysis toolbox RSLA! For a general intro to the toolbox, type "man rsla". But before jumping right into it, you should know how to use the manual. This will just take a minute!';

        choice=questdlg(['' WelcomeText ''],'WELCOME!','Manual Tutorial','Cancel','Manual Tutorial');

        switch choice
            case 'Manual Tutorial'
                msgbox('All interaction with the manual happens in the MATLAB console. Oh look!');
                disp('While a lot of information about the toolbox can be obtained by hovering over fields and buttons in the GUIs, the manual is supposed to give detailed information about how to use the toolbox.');
                fprintf('In some manual entries, you will see an arrow pointing to the right, directly followed by a word (e.g. ->rsla). This means that \nthat word has a manual entry which you can consult by typing "man TOPIC" (in this case: "man rsla").\n');
                fprintf('That''s it! In order to get started, just summon one of the GUIs you want to work with. If you don''t know how to do this, you should read the general information (->rsla).\n');
        end
    end

end