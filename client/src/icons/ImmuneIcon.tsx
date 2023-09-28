import Icon from './svg/immune_system.svg';

function ImmuneIcon({ size }: { size: number }) {
  return (
    <div style={{ width: size }}>
      <img src={Icon} alt="immune system icon" />
    </div>
  );
}

export default ImmuneIcon;
